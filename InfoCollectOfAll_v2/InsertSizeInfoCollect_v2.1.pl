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
my ($Bam,$File4InsertInfo,$MinSize,$MaxSize,$File4Cyto,$Bed,$File4Mark,$Mode4SizeInfo,$Flag4Single,$Flag4EndPos,$Flag4ListCheck,$Flag4Area,$ExtendLen,$Dir,$Samtools,$Bedtools);
my $HelpInfo = <<USAGE;

 $ThisScriptName
 Auther: zhangdong_xie\@foxmail.com

  This script can be used to collect insert size info.
  
  本脚本主要用于收集插入片段的长度信息。
  
  和v1版本相比，修复了一些bug。该版本可以当作正式版。
  v1版本本质上统计的是正向的片段长度信息（bam中第9列为正）。
  本版本会默认统计正/负的所有信息并去重，也可以指定统计正或者负的信息。
  
  和v2.0相比：
    1. 针对单点/区间统计新增一个模式，当选中时将由只统计长度数量，改为同时记录长度及对应起始坐标，便于回溯基因组序列信息。

 -i      ( Optional ) Bam file;
 -o      ( Optional ) File for result;

 # 长度范围限制
 -min    ( Optional ) Minimal size for count (100);
 -max    ( Optional ) Maximal size for count (250);

 # Bed文件相关
 -b      ( Optional ) File in bed format;
                      Bed文件，默认是3列。指定bed时按bed区域统计，否则统计全基因组。
 -snp    ( Optional ) If need to confirm the bed was snp list;
                      假如需要确认bed文件是否为对应的snp列表；
 -elen   ( Optional ) The extended length for snp list;
                      用于延申Bed区域的长度，针对收集那些应该覆盖了但没有测到的部分。
 
 # 按单点统计
 -s      ( Optional ) If collect info one point by one point (with -s only, only effective when -b specified);
                      假如需要对bed文件中的点逐个收集信息。
 
 # 按区间统计
 -a      ( Optional ) If collect info by area (with -a only, effective when -b specified);
                      假如按区间统计片段信息，不指定‘-m’时一行bed文件为统计单位，当指定‘-m’时以指定的第4列标签为统计单位。
 -m      ( Optional ) File of Mark for each line in Bed;
                      文件包括4列：染色体号、起始坐标、结束坐标、所属区域标记（一个标记可以对应多个bed区间）。
                      注意：该文件的前3列需要和bed文件保持一致。
 
 # 片段信息的统计维度，正、负或者都要
 -mode   ( Optional ) Mode for the size collection type, + or - or +/- (left, right or all ,default: all);
                      用于指定片段信息的统计模式，只统计正的（left）、负的（right）或者二者都有（all，默认）；
 -end    ( Optional ) If collect info for 5\'End position (with -pos only, only effective when -s/-a specified);
                      假如在统计片段长度的同时，还需要统计每条插入片段对应的起始坐标。指定‘-end’即可。
 
 -bin    ( Optional ) List for searching of related bin or scripts; 
 -h      ( Optional ) Help infomation;

USAGE

GetOptions(
	'i=s' => \$Bam,
	'o=s' => \$File4InsertInfo,
	'min:i' => \$MinSize,
	'max:i' => \$MaxSize,
	'b:s' => \$Bed,
	'snp!' => \$Flag4ListCheck,
	'elen:i' => \$ExtendLen,
	's!' => \$Flag4Single,
	'a!' => \$Flag4Area,
	'm:s' => \$File4Mark,
	'mode:s' => \$Mode4SizeInfo,
	'end!' => \$Flag4EndPos,
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
	$Mode4SizeInfo = "all" unless($Mode4SizeInfo);
	die "[ Error ] -mode not correct ($Mode4SizeInfo).\n" unless($Mode4SizeInfo eq "left" || $Mode4SizeInfo eq "right" || $Mode4SizeInfo eq "all");
	
	$BinList = BinListGet() if(!$BinList);
	$MinSize = BinSearch("MinSize",$BinList,1) unless($MinSize);
	$MaxSize = BinSearch("MaxSize",$BinList,1) unless($MaxSize);
	$Samtools = BinSearch("Samtools",$BinList);
	$File4Cyto = BinSearch("File4Cyto",$BinList);
}

if(1)
{
	# 检查是否为snplist文件，假如不是则用snplist文件替代；
	if($Flag4ListCheck && $Bed)
	{
		my $tBed = SNPListOnBed($Bed);
		die "[ Error ] Could not locate the snplist ($tBed) for $Bed.\n" unless($tBed && -s $tBed);
		$Bed = $tBed;
	}
	
	# 对bed延申、排序、合并以及调整区域定位文件顺序；
	my ($MergeBed,$File4Mark2) = ();
	# 不指定Bed文件时直接忽略
	if($Bed)
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
		printf "[ %s ] Bed file sort and merge done.\n",TimeString(time,$BeginTime);
		
		# 只有按区域统计时标签文件才起作用。
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
				die "[ Error ] Unknown mark for $Line\n" unless($Mark);
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
		my %SizeInfo = %{SizeInfoOfSingleOnBed2($Bam,$Samtools,$MergeBed,$MinSize,$MaxSize,$Mode4SizeInfo,$Flag4EndPos)};
		
		# 逐个点整理输出；
		open(BED,"cat $MergeBed |") or die $! unless($MergeBed =~ /\.gz$/);
		open(BED,"zcat $MergeBed |") or die $! if($MergeBed =~ /\.gz$/);
		open(IN,"> $File4InsertInfo") or die $! unless($File4InsertInfo =~ /\.gz$/);
		open(IN,"| gzip > $File4InsertInfo") or die $! if($File4InsertInfo =~ /\.gz$/);
		print IN join("\t","#Chr","Pos",@Size),"\n" unless($Flag4EndPos);
		print IN join("\t","#Chr","Pos","StringForSize","StringFor5Pos"),"\n" if($Flag4EndPos);
		while(my $Line = <BED>)
		{
			chomp $Line;
			my @Cols = split /\t/, $Line;
			for my $i ($Cols[1] + 1 .. $Cols[2])
			{
				my $Key = join("\t",$Cols[0],$i);
				
				if($Flag4EndPos)
				{
					my (@SSize,@SPos) = ();
					for my $j (0 .. $#{$SizeInfo{$Key}})
					{
						my @SItems = split /,/, $SizeInfo{$Key}[$j];
						push @SSize, $SItems[1];
						push @SPos, $SItems[0];
					}
					my $String4Size = join(",",@SSize);
					my $String4Pos = join(",",@SPos);
					print IN join("\t",$Key,$String4Size,$String4Pos),"\n";
				}
				else
				{
					my @SingleSizeInfo = ();
					for my $j (0 .. $SizeRange)
					{
						$SizeInfo{$Key}[$j] = 0 unless($SizeInfo{$Key}[$j]);
						$SingleSizeInfo[$j] = $SizeInfo{$Key}[$j];
					}
					print IN join("\t",$Key,@SingleSizeInfo),"\n";
				}
			}
		}
		close BED;
		close IN;
	}
	# 对Bed文件按区域统计；
	elsif($MergeBed && $Flag4Area)
	{
		# 假如不指定Mark分类标签，则按Bed行统计；
		my $Flag4Mark = "";
		$Flag4Mark = $File4Mark2 if($File4Mark2);
		my %SizeInfo = %{SizeInfoOfMultiArea2($Bam,$Samtools,$MergeBed,$Flag4Mark,$MinSize,$MaxSize,$Mode4SizeInfo,$Flag4EndPos)};
		open(IN,"> $File4InsertInfo") or die $! unless($File4InsertInfo =~ /\.gz$/);
		open(IN,"| gzip > $File4InsertInfo") or die $! if($File4InsertInfo =~ /\.gz$/);
		if($Flag4Mark)
		{
			print IN join("\t","#Mark",@Size),"\n" unless($Flag4EndPos);
			print IN join("\t","#Mark","StringForSize","StringFor5Pos"),"\n" if($Flag4EndPos);
		}
		else
		{
			print IN join("\t","#Chr","From","To",@Size),"\n" unless($Flag4EndPos);
			print IN join("\t","#Chr","From","To","StringForSize","StringFor5Pos"),"\n" if($Flag4EndPos);
		}
		foreach my $Mark (keys %SizeInfo)
		{
			my $CMark = $Mark;
			unless($Flag4Mark)
			{
				my @Items = split /_/, $Mark;
				$CMark = join("\t",@Items);
			}
			
			if($Flag4EndPos)
			{
				my (@SSize,@SPos) = ();
				for my $j (0 .. $#{$SizeInfo{$Mark}})
				{
					my @SItems = split /,/, $SizeInfo{$Mark}[$j];
					push @SSize, $SItems[1];
					push @SPos, $SItems[0];
				}
				my $String4Size = join(",",@SSize);
				my $String4Pos = join(",",@SPos);
				print IN join("\t",$CMark,$String4Size,$String4Pos),"\n";
			}
			else
			{
				print IN join("\t",$CMark,@{$SizeInfo{$Mark}}),"\n";
			}
		}
		close IN;
		`rm $Flag4Mark` if($Flag4Mark && -s $Flag4Mark);
	}
	# 或者按染色体分析，指定或不指定bed文件均可;
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
