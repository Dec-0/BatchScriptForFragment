#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
Getopt::Long::Configure qw(no_ignore_case);
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin/.Modules";
use Parameter::BinList;
use MGDRelated::InsertSize;
use MGDRelated::InsertSize2;
use Sort::ChrPos;
use BedRelated::Bed;

my ($HelpFlag,$BinList,$BeginTime);
my $ThisScriptName = basename $0;
my ($Bam,$File4InsertInfo,$MinSize,$MaxSize,$File4Cyto,$Bed,$File4Mark,$Flag4Single,$Flag4ListCheck,$Flag4Area,$ExtendLen,$Dir,$Samtools,$Bedtools);
my $HelpInfo = <<USAGE;

 $ThisScriptName
 Auther: zhangdong_xie\@foxmail.com

  This script was used to collect insert size info.
  本脚本主要用于收集插入片段的长度信息。
  
  和v1版本相比，修改如下：
    1. 重新修改了单点及区域数据收集函数。
    2. 新增区域标记选项。

 -i      ( Optional ) Bam file;
 -o      ( Optional ) File for result;

 -min    ( Optional ) Minimal size for count (100);
 -max    ( Optional ) Maximal size for count (250);
 -b      ( Optional ) Bed;
 -m      ( Optional ) File of Mark for each line in Bed;
                      文件包括4列：染色体号、起始坐标、结束坐标、所属区域标记（一个标记可以对应多个bed区间）。
                      注意：该文件的前3列需要和bed文件保持一致。
 -s      ( Optional ) If collect info one point by one point (with -s only, only effective when -b specified);
                      假如需要对bed文件中的点逐个收集信息。
 -snp    ( Optional ) If need to confirm the bed was snp list;
                      假如需要确认bed文件是否为对应的snp列表；
 -a      ( Optional ) If collect info by area (with -a only, effective when -b specified);
                      假如按区间统计；
 -elen   ( Optional ) The extended length for snp list;
                      用于延申Bed区域的长度，针对收集那些应该覆盖了但没有测到的部分。
 -bin    ( Optional ) List for searching of related bin or scripts; 
 -h      ( Optional ) Help infomation;

USAGE

GetOptions(
	'i=s' => \$Bam,
	'o=s' => \$File4InsertInfo,
	'min:i' => \$MinSize,
	'max:i' => \$MaxSize,
	'b:s' => \$Bed,
	'm:s' => \$File4Mark,
	's!' => \$Flag4Single,
	'snp!' => \$Flag4ListCheck,
	'a!' => \$Flag4Area,
	'elen:i' => \$ExtendLen,
	'bin:s' => \$BinList,
	'h!' => \$HelpFlag
) or die $HelpInfo;

if($HelpFlag || !$Bam || !$File4InsertInfo)
{
	die $HelpInfo;
}
else
{
	$BeginTime = ScriptBegin(0,$ThisScriptName);
	IfFileExist($Bam);
	IfFileExist($Bed) if($Bed);
	$Bed = "" unless($Bed);
	$Dir = dirname $File4InsertInfo;
	IfDirExist($Dir);
	$ExtendLen = 0 unless($ExtendLen);
	
	$BinList = BinListGet() if(!$BinList);
	$Samtools = BinSearch("Samtools",$BinList);
	$MinSize = BinSearch("MinSize",$BinList,1);
	$MaxSize = BinSearch("MaxSize",$BinList,1);
	$File4Cyto = BinSearch("File4Cyto",$BinList,1);
}

if(1)
{
	# 检查是否为snplist文件，假如不是则用snplist文件替代；
	if($Flag4ListCheck && $Bed)
	{
		my $tBed = SNPListOnBed($Bed);
		die "[ Error ] Could not locate the snplist for $Bed.\n" unless($tBed && -s $tBed);
		$Bed = $tBed;
	}
	
	# 对bed延申、排序、合并以及调整区域定位文件顺序；
	my ($MergeBed,$File4Mark2) = ();
	if(1)
	{
		$Bedtools = BinSearch("Bedtools",$BinList);
		
		# 对Bed进行排序；
		my (@Chr,@From,@To,@OtherInfo) = ();
		open(BED,"cat $Bed |") or die $! unless($Bed =~ /\.gz$/);
		open(BED,"zcat $Bed |") or die $! if($Bed =~ /\.gz$/);
		while(my $Line = <BED>)
		{
			chomp $Line;
			my @Cols = split /\t/, $Line;
			if($ExtendLen)
			{
				$Cols[1] -= $ExtendLen;
				$Cols[2] += $ExtendLen;
			}
			push @Chr, $Cols[0];
			push @From, $Cols[1];
			push @To, $Cols[2];
			push @OtherInfo, "";
		}
		close BED;
		my @Ref = ChrPosAndOtherByNum(\@Chr,\@From,\@To,\@OtherInfo);
		@Chr = @{$Ref[0]};
		@From = @{$Ref[1]};
		@To = @{$Ref[2]};
		my $BaseName = basename $Bed;
		$BaseName =~ s/\.gz$//;
		$BaseName =~ s/\.bed$//;
		my $RandomString = RandString(10);
		my $SortBed = $Dir . "/" . $BaseName . ".Sort." . $RandomString . ".bed";
		open(SORT,"> $SortBed") or die $!;
		for my $i (0 .. $#Chr)
		{
			print SORT join("\t",$Chr[$i],$From[$i],$To[$i]),"\n";
		}
		close SORT;
		
		# 对区域进行合并并提取数据；
		$MergeBed = $Dir . "/" . $BaseName . ".Sort.Merge." . $RandomString . ".bed";
		`$Bedtools merge -i $SortBed > $MergeBed`;
		`rm $SortBed` if($SortBed && -s $SortBed);
		
		if($File4Mark && $Flag4Area)
		{
			$File4Mark2 = $Dir . "/" . $BaseName . ".Sort.Merge." . $RandomString . ".mark";
			open(MERGE,"cat $MergeBed |") or die $!;
			open(MARK,"> $File4Mark2") or die $!;
			while(my $Line = <MERGE>)
			{
				chomp $Line;
				my @Cols = split /\t/, $Line;
				my $Mark = "";
				$Mark = `cat $File4Mark | awk '{if(\$1 == "$Cols[0]" && \$2 <= $Cols[2] && \$3 >= $Cols[1]){print \$4}}' | head -n1` unless($File4Mark =~ /\.gz$/);
				$Mark = `zcat $File4Mark | awk '{if(\$1 == "$Cols[0]" && \$2 <= $Cols[2] && \$3 >= $Cols[1]){print \$4}}' | head -n1` if($File4Mark =~ /\.gz$/);
				chomp $Mark;
				print MARK join("\t",@Cols[0 .. 2],$Mark),"\n";
			}
			close MERGE;
			close MARK;
		}
	}
	
	my @Size = ($MinSize .. $MaxSize);
	my $SizeRange = $MaxSize - $MinSize;
	# 对Bed文件，按点统计；
	if($MergeBed && $Flag4Single)
	{
		# 统计所有信息；
		my %SizeInfo = %{SizeInfoOfSingleOnBed2($Bam,$Samtools,$MergeBed,$MinSize,$MaxSize)};
		
		# 逐个点整理输出；
		open(BED,"cat $MergeBed |") or die $! unless($MergeBed =~ /\.gz$/);
		open(BED,"zcat $MergeBed |") or die $! if($MergeBed =~ /\.gz$/);
		open(IN,"> $File4InsertInfo") or die $! unless($File4InsertInfo =~ /\.gz$/);
		open(IN,"| gzip > $File4InsertInfo") or die $! if($File4InsertInfo =~ /\.gz$/);
		print IN join("\t","#Chr","Pos",@Size),"\n";
		while(my $Line = <BED>)
		{
			chomp $Line;
			my @Cols = split /\t/, $Line;
			for my $i ($Cols[1] + 1 .. $Cols[2])
			{
				my $Key = join("\t",$Cols[0],$i);
				my @SingleSizeInfo = ();
				for my $j (0 .. $SizeRange)
				{
					$SizeInfo{$Key}[$j] = 0 unless($SizeInfo{$Key}[$j]);
					$SingleSizeInfo[$j] = $SizeInfo{$Key}[$j];
				}
				
				print IN join("\t",$Key,@SingleSizeInfo),"\n";
			}
		}
		close BED;
		close IN;
	}
	# 对Bed文件按区域统计；
	elsif($MergeBed && $Flag4Area)
	{
		my $Flag4Mark = "";
		$Flag4Mark = $File4Mark2 if($File4Mark2);
		my %SizeInfo = %{SizeInfoOfMultiArea2($Bam,$Samtools,$MergeBed,$Flag4Mark,$MinSize,$MaxSize)};
		open(IN,"> $File4InsertInfo") or die $! unless($File4InsertInfo =~ /\.gz$/);
		open(IN,"| gzip > $File4InsertInfo") or die $! if($File4InsertInfo =~ /\.gz$/);
		print IN join("\t","#Mark",@Size),"\n";
		foreach my $Mark (keys %SizeInfo)
		{
			print IN join("\t",$Mark,@{$SizeInfo{$Mark}}),"\n";
		}
		close IN;
		`rm $Flag4Mark` if($Flag4Mark && -s $Flag4Mark);
	}
	# 没有Bed文件时，按染色体分析;
	else
	{
		my %InsertInfo = %{SizeInfoGetOnChr($Bam,$Samtools,$MergeBed,$MinSize,$MaxSize)};
		open(IN,"> $File4InsertInfo") or die $! unless($File4InsertInfo =~ /\.gz$/);
		open(IN,"| gzip > $File4InsertInfo") or die $! if($File4InsertInfo =~ /\.gz$/);
		print IN join("\t","#Chr",@Size),"\n";
		foreach my $Chr (keys %InsertInfo)
		{
			print IN join("\t",$Chr,@{$InsertInfo{$Chr}}),"\n";
		}
		close IN;
	}
	
	`rm $MergeBed` if($MergeBed && -s $MergeBed);
}
printf "[ %s ] The end.\n",TimeString(time,$BeginTime);


######### Sub functions ##########
