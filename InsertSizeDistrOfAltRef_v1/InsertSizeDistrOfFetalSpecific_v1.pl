#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
Getopt::Long::Configure qw(no_ignore_case);
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin/.Modules";
use Parameter::BinList;

my ($HelpFlag,$BinList,$BeginTime);
my $ThisScriptName = basename $0;
my ($BedWithMark,$File4GTypePair,$File4RefIS,$File4AltIS,$File4SizeInfo,$MinSize,$MaxSize,$Mode4EndPrefer);
my @String4PairInfo;
my $HelpInfo = <<USAGE;

 $ThisScriptName
 Auther: zhangdong_xie\@foxmail.com

  This script can be used to count the size number of specific maternal-fetal geno-type pair in bed area.
  
  本脚本主要用于统计胎儿/母本特异部分的片段数量信息.
  假如统计胎儿特异的部分可以指定“AA,AB,Ref BB,BA,Alt”，统计母本特异的可以指定“AB,AA,Ref AB,BB,Alt”（其中A代表Alt，B代表Ref）。

 -g      ( Required ) File with at least 7 columns "chr,from,to,ref,alt,maternalGType,fetalGType";
                      记录位点和母子基因型的文件。
 -ref    ( Required ) File for the insert info counts of reference;
                      Ref碱基对应的特定插入片段信息。文件至少有7列：染色体号、起始坐标、结束坐标、ref序列、当前序列（ref或alt）、片段长度信息（逗号分隔）、位点相对5\'位置信息。
 -alt    ( Required ) File for the insert info counts of allele;
                      Alt碱基对应的特定插入片段信息；
 -o      ( Required ) File for result;

 # 统计区域及基因型分类；
 -b      ( Required ) Bed with mark in the 4th column;
                      限制统计的范围以及区域的分类标签。bed文件第4列为分类标签。
 -p      ( Required ) PairInfo for collecting, with format like 'AA,AB,Ref', 'BB,BA,Alt' (Multi times);
                      指定需要统计的母子配对基因型，以及该配对基因型需要纳入统计的Ref或Alt标签。

 # 插入片段长度统计范围
 -min    ( Optional ) Minimal size for count (50);
 -max    ( Optional ) Maximal size for count (300);
 
 # 片段信息的筛选参数
 -mode   ( Optional ) Mode for which end the site prefer (left, right or both, default: both);
                      变异位置与片段末端的倾向性筛选。SNP位点靠近片段左侧、右侧或者都行。

 -bin    ( Optional ) List for searching of related bin or scripts; 
 -h      ( Optional ) Help infomation;

USAGE

GetOptions(
	'b=s' => \$BedWithMark,
	'g=s' => \$File4GTypePair,
	'ref=s' => \$File4RefIS,
	'alt=s' => \$File4AltIS,
	'p=s' => \@String4PairInfo,
	'o=s' => \$File4SizeInfo,
	'min:i' => \$MinSize,
	'max:i' => \$MaxSize,
	'mode:s' => \$Mode4EndPrefer,
	'bin:s' => \$BinList,
	'h!' => \$HelpFlag
) or die $HelpInfo;

if($HelpFlag || !$BedWithMark || !$File4GTypePair || !$File4RefIS || !$File4AltIS || !@String4PairInfo)
{
	die $HelpInfo;
}
else
{
	$BeginTime = ScriptBegin(0,$ThisScriptName);
	IfFileExist($BedWithMark,$File4GTypePair,$File4RefIS,$File4AltIS);
	
	$BinList = BinListGet() if(!$BinList);
	$MinSize = BinSearch("MinSize",$BinList,1) unless($MinSize);
	$MaxSize = BinSearch("MaxSize",$BinList,1) unless($MaxSize);
	die "[ Error ] MinSize > MaxSize ($MinSize > $MaxSize)\n" if($MinSize > $MaxSize);
	$Mode4EndPrefer = "both" unless($Mode4EndPrefer);
	die "[ Error ] Mode for End Prefer not correct ($Mode4EndPrefer)\n" unless($Mode4EndPrefer eq "left" || $Mode4EndPrefer eq "right" || $Mode4EndPrefer eq "both");
}

if(1)
{
	my (%PosInfo,%Mark) = ();
	if(1)
	{
		# 记录位点对应的区域标记；
		my %PosMark = ();
		open(BM,"cat $BedWithMark | grep -v ^# |") or die $! unless($BedWithMark =~ /\.gz$/);
		open(BM,"zcat $BedWithMark | grep -v ^# |") or die $! if($BedWithMark =~ /\.gz$/);
		while(my $Line = <BM>)
		{
			chomp $Line;
			my @Cols = split /\t/, $Line;
			for my $i ($Cols[1] + 1 .. $Cols[2])
			{
				my $Key = join("\t",$Cols[0],$i);
				$PosMark{$Key} = $Cols[3];
				$Mark{$Cols[3]} = 1;
			}
		}
		close BM;
		
		# 记录母子配对基因型对应的文件类型，比如Ref或Alt；
		my %FileMark = ();
		for my $i (0 .. $#String4PairInfo)
		{
			my ($MGType,$FGType,$tFileMark) = split /,/, $String4PairInfo[$i];
			$FileMark{$MGType}{$FGType} = $tFileMark;
			die "[ Error ] File Mark not Ref or Alt ($String4PairInfo[$i])\n" unless($tFileMark eq "Ref" || $tFileMark eq "Alt");
		}
		
		# 记录位点对应的Ref或Alt属性，和区域标记的关系；
		open(GTP,"cat $File4GTypePair | grep -v ^# |") or die $! unless($File4GTypePair =~ /\.gz$/);
		open(GTP,"zcat $File4GTypePair | grep -v ^# |") or die $! if($File4GTypePair =~ /\.gz$/);
		while(my $Line = <GTP>)
		{
			chomp $Line;
			my @Cols = split /\t/, $Line;
			# 排除掉异常基因型对应的位点
			next unless($FileMark{$Cols[5]}{$Cols[6]});
			my $Key = join("\t",@Cols[0 .. 1]);
			# 排除掉Bed区域外的点
			next unless($PosMark{$Key});
			$PosInfo{$FileMark{$Cols[5]}{$Cols[6]}}{$Key} = $PosMark{$Key};
		}
		close GTP;
	}
	
	my %SizeInfo = ();
	if(1)
	{
		open(REF,"cat $File4RefIS | grep -v ^# |") or die $! unless($File4RefIS =~ /\.gz$/);
		open(REF,"zcat $File4RefIS | grep -v ^# |") or die $! if($File4RefIS =~ /\.gz$/);
		while(my $Line = <REF>)
		{
			chomp $Line;
			my @Cols = split /\t/, $Line;
			next if($Cols[5] eq "-");
			my $Key = join("\t",@Cols[0 .. 1]);
			next unless($PosInfo{"Ref"}{$Key});
			
			my @Size = split /,/, $Cols[5];
			my @Pos = split /,/, $Cols[6];
			for my $i (0 .. $#Size)
			{
				next if($Size[$i] < $MinSize || $Size[$i] > $MaxSize);
				if($Mode4EndPrefer eq "left")
				{
					next if($Pos[$i] * 2 > $Size[$i] && $Pos[$i] > 100);
				}
				elsif($Mode4EndPrefer eq "right")
				{
					next if($Pos[$i] * 2 < $Size[$i] && $Size[$i] > $Pos[$i] + 100);
				}
				my $tId = $Size[$i] - $MinSize;
				$SizeInfo{$PosInfo{"Ref"}{$Key}}[$tId] = 0 unless($SizeInfo{$PosInfo{"Ref"}{$Key}}[$tId]);
				$SizeInfo{$PosInfo{"Ref"}{$Key}}[$tId] ++;
			}
		}
		close REF;
		
		
		open(ALT,"cat $File4AltIS | grep -v ^# |") or die $! unless($File4AltIS =~ /\.gz$/);
		open(ALT,"zcat $File4AltIS | grep -v ^# |") or die $! if($File4AltIS =~ /\.gz$/);
		while(my $Line = <ALT>)
		{
			chomp $Line;
			my @Cols = split /\t/, $Line;
			next if($Cols[5] eq "-");
			my $Key = join("\t",@Cols[0 .. 1]);
			next unless($PosInfo{"Alt"}{$Key});
			
			my @Size = split /,/, $Cols[5];
			my @Pos = split /,/, $Cols[6];
			for my $i (0 .. $#Size)
			{
				next if($Size[$i] < $MinSize || $Size[$i] > $MaxSize);
				if($Mode4EndPrefer eq "left")
				{
					next if($Pos[$i] * 2 > $Size[$i] && $Pos[$i] > 100);
				}
				elsif($Mode4EndPrefer eq "right")
				{
					next if($Pos[$i] * 2 < $Size[$i] && $Size[$i] > $Pos[$i] + 100);
				}
				my $tId = $Size[$i] - $MinSize;
				$SizeInfo{$PosInfo{"Alt"}{$Key}}[$tId] = 0 unless($SizeInfo{$PosInfo{"Alt"}{$Key}}[$tId]);
				$SizeInfo{$PosInfo{"Alt"}{$Key}}[$tId] ++;
			}
		}
		close ALT;
	}
	
	open(SI,"> $File4SizeInfo") or die $! unless($File4SizeInfo =~ /\.gz$/);
	open(SI,"| gzip > $File4SizeInfo") or die $! if($File4SizeInfo =~ /\.gz$/);
	my @Size = ($MinSize .. $MaxSize);
	print SI join("\t","#Mark",@Size),"\n";
	foreach my $Key (keys %Mark)
	{
		for my $i (0 .. $#Size)
		{
			$SizeInfo{$Key}[$i] = 0 unless($SizeInfo{$Key}[$i]);
		}
		print SI join("\t",$Key,@{$SizeInfo{$Key}}),"\n";
	}
	close SI;
}
printf "[ %s ] The end.\n",TimeString(time,$BeginTime);


######### Sub functions ##########
