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
use Sort::ChrPos;
use BedRelated::Bed;

my ($HelpFlag,$BinList,$BeginTime);
my $ThisScriptName = basename $0;
my ($Bam,$File4InsertInfo,$MinSize,$MaxSize,$Bed,$Flag4Single,$Flag4ListCheck,$Flag4Area,$ExtendLen,$Dir,$Samtools,$Bedtools);
my $HelpInfo = <<USAGE;

 $ThisScriptName
 Auther: zhangdong_xie\@foxmail.com

  This script was used to collect insert size info.
  本脚本主要用于收集插入片段的长度信息。

 -i      ( Optional ) Bam file;
 -o      ( Optional ) File for result;

 -min    ( Optional ) Minimal size for count (100);
 -max    ( Optional ) Maximal size for count (250);
 -b      ( Optional ) Bed;
 -s      ( Optional ) If collect info one point by one point (with -s only, only effective when -b specified);
                      假如需要对bed文件中的点逐个收集信息。
 -snp    ( Optional ) If need to confirm the bed was snp list;
                      假如需要确认bed文件是否为对应的snp列表；
 -a      ( Optional ) If collect info by area (with -a only, effective when -b specified);
                      假如按区间统计；
 -elen   ( Optional ) The extended length for snp list;
                      用于指定snplist列表延申的长度；
 -bin    ( Optional ) List for searching of related bin or scripts; 
 -h      ( Optional ) Help infomation;

USAGE

GetOptions(
	'i=s' => \$Bam,
	'o=s' => \$File4InsertInfo,
	'min:i' => \$MinSize,
	'max:i' => \$MaxSize,
	'b:s' => \$Bed,
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
}

if(1)
{
	# 检查是否为snplist文件；
	if($Flag4ListCheck && $Bed)
	{
		my $tBed = SNPListOnBed($Bed);
		die "[ Error ] Could not locate the snplist for $Bed.\n" unless($tBed && -s $tBed);
		$Bed = $tBed;
	}
	
	my @Size = ($MinSize .. $MaxSize);
	my $SizeRange = $MaxSize - $MinSize;
	# 统计Bed文件内每个点对应的片段长度分布信息；
	if($Bed && $Flag4Single)
	{
		# 统计所有信息；
		my %SizeInfo = %{SizeInfoOfSingleOnBed($Bam,$Samtools,$Bed,$MinSize,$MaxSize)};
		
		# 逐个点整理输出；
		open(BED,"cat $Bed |") or die $! unless($Bed =~ /\.gz$/);
		open(BED,"zcat $Bed |") or die $! if($Bed =~ /\.gz$/);
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
	# 对Bed文件每行，逐行统计；
	elsif($Bed && $Flag4Area)
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
		my $MergeBed = $Dir . "/" . $BaseName . ".Sort.Merge." . $RandomString . ".bed";
		`$Bedtools merge -i $SortBed > $MergeBed`;
		my @SizeInfo = @{SizeInfoOfMultiArea($Bam,$Samtools,$SortBed,$MergeBed,$MinSize,$MaxSize)};
		open(IN,"> $File4InsertInfo") or die $! unless($File4InsertInfo =~ /\.gz$/);
		open(IN,"| gzip > $File4InsertInfo") or die $! if($File4InsertInfo =~ /\.gz$/);
		print IN join("\t","#Chr","From","To",@Size),"\n";
		for my $i (0 .. $#Chr)
		{
			print IN join("\t",$Chr[$i],$From[$i],$To[$i],@{$SizeInfo[$i]}),"\n";
		}
		close IN;
		`rm $SortBed` if(-s $SortBed);
		`rm $MergeBed` if(-s $MergeBed);
	}
	# 假如不给Bed文件就整体分析;
	else
	{
		my %InsertInfo = %{SizeInfoGetOnChr($Bam,$Samtools,$Bed,$MinSize,$MaxSize)};
		open(IN,"> $File4InsertInfo") or die $! unless($File4InsertInfo =~ /\.gz$/);
		open(IN,"| gzip > $File4InsertInfo") or die $! if($File4InsertInfo =~ /\.gz$/);
		print IN join("\t","#Chr",@Size),"\n";
		foreach my $Chr (keys %InsertInfo)
		{
			print IN join("\t",$Chr,@{$InsertInfo{$Chr}}),"\n";
		}
		close IN;
	}
}
printf "[ %s ] The end.\n",TimeString(time,$BeginTime);


######### Sub functions ##########
