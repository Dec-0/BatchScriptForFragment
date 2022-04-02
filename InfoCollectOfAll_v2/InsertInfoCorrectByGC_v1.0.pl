#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
Getopt::Long::Configure qw(no_ignore_case);
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin/.Modules";
use Parameter::BinList;
use SeqRelated::Seq;

my ($HelpFlag,$BinList,$BeginTime);
my $ThisScriptName = basename $0;
my ($File4OriInfo,$File4CorrectedInfo,$MinSize,$MaxSize,$Reference,$Bedtools);
my $HelpInfo = <<USAGE;

 $ThisScriptName
 Auther: zhangdong_xie\@foxmail.com

  This script can be used to correct the number of fragments by gc info of reference.
  
  本脚本用于处理InsertSizeInfoCollect在‘-end’模式下生成的文件，利用GC信息对片段分布进行校准。
  采用的策略是Benjamini et al., 2012介绍的fragment model。

 -i      ( Required ) File generated by InsertSizeInfoCollect, with \"StringForSize\" and \"StringFor5Pos\";
                      输入文件前三列为探针信息，后两列分别为frgment长度以及对应的5’起始位置。
 -o      ( Required ) File for the corrected result;

 -min    ( Optional ) Minimal size for count (100);
 -max    ( Optional ) Maximal size for count (250);
 
 -bin    ( Optional ) List for searching of related bin or scripts; 
 -h      ( Optional ) Help infomation;

USAGE

GetOptions(
	'i=s' => \$File4OriInfo,
	'o=s' => \$File4CorrectedInfo,
	'min:i' => \$MinSize,
	'max:i' => \$MaxSize,
	'bin:s' => \$BinList,
	'h!' => \$HelpFlag
) or die $HelpInfo;

if($HelpFlag || !$File4OriInfo || !$File4CorrectedInfo)
{
	die $HelpInfo;
}
else
{
	$BeginTime = ScriptBegin(0,$ThisScriptName);
	
	
	$BinList = BinListGet() if(!$BinList);
	$MinSize = BinSearch("MinSize",$BinList,1) unless($MinSize);
	$MaxSize = BinSearch("MaxSize",$BinList,1) unless($MaxSize);
	$Reference = BinSearch("Reference",$BinList);
	$Bedtools = BinSearch("Bedtools",$BinList);
}

if(1)
{
	# 1. 过滤掉5'不在芯片区域内的以及长度不达标的数据；
	my ($File4Flt,$File4SizeStatBeforeCorrect) = ();
	my @SizeMark = ($MinSize .. $MaxSize);
	my $SizeNumTotal = 0;
	if($File4OriInfo)
	{
		$File4Flt = $File4CorrectedInfo;
		$File4Flt =~ s/\.gz$//;
		$File4Flt =~ s/\.xls$//;
		$File4SizeStatBeforeCorrect = $File4Flt;
		$File4Flt .= ".PosFlt.xls.gz";
		$File4SizeStatBeforeCorrect .= ".BeforeCorrect.xls.gz";
		
		my $SizeRangeLen = $MaxSize - $MinSize + 1;
		open(ORI,"cat $File4OriInfo |") or die $! unless($File4OriInfo =~ /\.gz$/);
		open(ORI,"zcat $File4OriInfo |") or die $! if($File4OriInfo =~ /\.gz$/);
		open(FLT,"> $File4Flt") or die $! unless($File4Flt =~ /\.gz$/);
		open(FLT,"| gzip > $File4Flt") or die $! if($File4Flt =~ /\.gz$/);
		open(STAT,"> $File4SizeStatBeforeCorrect") or die $! unless($File4SizeStatBeforeCorrect =~ /\.gz$/);
		open(STAT,"| gzip > $File4SizeStatBeforeCorrect") or die $! if($File4SizeStatBeforeCorrect =~ /\.gz$/);
		while(my $Line = <ORI>)
		{
			if($Line =~ /^#/)
			{
				print FLT $Line;
				print STAT join("\t","#Chr","From","To",@SizeMark),"\n";
				next;
			}
			chomp $Line;
			my ($Chr,$From,$To,$String4Size,$String4Pos) = split /\t/, $Line;
			my @Size = split /,/, $String4Size;
			my @Pos = split /,/, $String4Pos;
			
			my (@tSize,@tPos) = ();
			for my $i (0 .. $#Pos)
			{
				next if($Pos[$i] < $From || $Pos[$i] > $To);
				next if($Size[$i] < $MinSize || $Size[$i] > $MaxSize);
				push @tSize, $Size[$i];
				push @tPos, $Pos[$i];
			}
			next unless(@tSize);
			$String4Size = join(",",@tSize);
			$String4Pos = join(",",@tPos);
			print FLT join("\t",$Chr,$From,$To,$String4Size,$String4Pos),"\n";
			
			my @SizeStat = (0) x $SizeRangeLen;
			for my $i (0 .. $#tSize)
			{
				$SizeStat[$tSize[$i] - $MinSize] ++;
				$SizeNumTotal ++;
			}
			print STAT join("\t",$Chr,$From,$To,@SizeStat),"\n";
		}
		close ORI;
		close FLT;
		close STAT;
	}
	printf "[ %s ] Done for flt.\n",TimeString(time,$BeginTime);
	
	# 2. 逐个点统计长度和gc组合对应的数量信息；
	my %OriValueOfSizePlusGC = ();
	if($File4Flt)
	{
		open(FLT,"cat $File4Flt |") or die $! unless($File4Flt =~ /\.gz$/);
		open(FLT,"zcat $File4Flt |") or die $! if($File4Flt =~ /\.gz$/);
		while(my $Line = <FLT>)
		{
			next if($Line =~ /^#/);
			chomp $Line;
			my ($Chr,$From,$To,$String4Size,$String4Pos) = split /\t/, $Line;
			my @Size = split /,/, $String4Size;
			my @Pos = split /,/, $String4Pos;
			
			my ($MinPos,$MaxPos) = ($Pos[0],$Pos[0]);
			for my $i (0 .. $#Pos)
			{
				$MinPos = $Pos[$i] if($Pos[$i] < $MinPos);
				my $tPos = $Pos[$i] + $Size[$i] - 1;
				$MaxPos = $tPos if($tPos > $MaxPos);
			}
			my $LSeq = RefGet($Reference,$Chr,$MinPos,$MaxPos,$Bedtools);
			
			# Wa,l 计算的是5’起始位置相同的片段数量，需要将坐标相同的合并；
			my %SizePosNum = ();
			for my $i (0 .. $#Size)
			{
				$SizePosNum{$Pos[$i]}{$Size[$i]} = 0 unless($SizePosNum{$Pos[$i]}{$Size[$i]});
				$SizePosNum{$Pos[$i]}{$Size[$i]} ++;
			}
			
			foreach my $KPos (keys %SizePosNum)
			{
				foreach my $KSize (keys %{$SizePosNum{$KPos}})
				{
					my $tFrom = $KPos - $MinPos + 2;
					my $tSize = $KSize - 4;
					my $Seq = substr($LSeq,$tFrom,$tSize);
					my $GCNum = GCNum($Seq);
					my $Key = join("\t",$KSize,$GCNum);
					push @{$OriValueOfSizePlusGC{$Key}}, $SizePosNum{$KPos}{$KSize};
				}
			}
		}
		close FLT;
	}
	printf "[ %s ] Done for size + gc combination.\n",TimeString(time,$BeginTime);
	
	# 3. 计算理论值以及记录；
	my %PValueOfSizePlusGC = ();
	if(%OriValueOfSizePlusGC)
	{
		my $File4Arg = $File4Flt;
		$File4Arg =~ s/\.PosFlt\.xls\.gz$//;
		$File4Arg .= ".FragmentModel.arg";
		open(LOG,"> $File4Arg") or die $!;
		print LOG join("\t","#Size","GCNum","PredictedValue","OriginalValue"),"\n";
		foreach my $Key (sort {$a cmp $b} keys %OriValueOfSizePlusGC)
		{
			my ($Size,$GCNum) = split /\t/, $Key;
			my @Num = sort {$a <=> $b} @{$OriValueOfSizePlusGC{$Key}};
			my $LastId = @Num;
			$LastId = int($LastId * 0.99 + 0.5);
			my @NNum = @Num[0 .. $LastId - 1];
			my $Mean = 0;
			for my $i (0 .. $#NNum)
			{
				$Mean += $NNum[$i];
			}
			$Mean = sprintf("%.8f",$Mean / $LastId);
			my $String4Num = join(",",@Num);
			print LOG join("\t",$Size,$GCNum,$Mean,$String4Num),"\n";
			$PValueOfSizePlusGC{$Key} = $Mean;
		}
		close LOG;
	}
	printf "[ %s ] Done for prediction value.\n",TimeString(time,$BeginTime);
	
	# 4. 逐个探针矫正
	my @SizeAll = ();
	if(%PValueOfSizePlusGC)
	{
		my $tId = 0;
		my $SizeRangeLen = $MaxSize - $MinSize + 1;
		open(FLT,"cat $File4Flt |") or die $! unless($File4Flt =~ /\.gz$/);
		open(FLT,"zcat $File4Flt |") or die $! if($File4Flt =~ /\.gz$/);
		while(my $Line = <FLT>)
		{
			next if($Line =~ /^#/);
			chomp $Line;
			my ($Chr,$From,$To,$String4Size,$String4Pos) = split /\t/, $Line;
			my @Size = split /,/, $String4Size;
			my @Pos = split /,/, $String4Pos;
			
			my $MaxPos = $To + $MaxSize;
			my $LSeq = RefGet($Reference,$Chr,$From,$MaxPos,$Bedtools);
			
			# 将序列信息转换成GC信息；
			my (@LBase,@LBase4GC) = ();
			my $AccumGCNum = 0;
			@LBase = split //, $LSeq;
			for my $i (0 .. $#LBase)
			{
				$AccumGCNum ++ if($LBase[$i] eq "G" || $LBase[$i] eq "C");
				$LBase4GC[$i] = $AccumGCNum;
			}
			
			@{$SizeAll[$tId]} = (0) x $SizeRangeLen;
			for my $i ($From .. $To)
			{
				for my $j ($MinSize .. $MaxSize)
				{
					my $tFrom = $i + 2 - $From;
					my $tTo = $tFrom + $j - 5;
					my $GCNum = $LBase4GC[$tTo] - $LBase4GC[$tFrom];
					
					my $Key = join("\t",$j,$GCNum);
					next unless($PValueOfSizePlusGC{$Key});
					$SizeAll[$tId][$i - $From] += $PValueOfSizePlusGC{$Key};
				}
			}
			$tId ++;
		}
		close FLT;
	}
	printf "[ %s ] Done for correction of single probe.\n",TimeString(time,$BeginTime);
	
	# 5. 总数校准；
	if($SizeNumTotal)
	{
		my $SizeNumTotalNew = 0;
		for my $i (0 .. $#SizeAll)
		{
			for my $j (0 .. $#{$SizeAll[$i]})
			{
				$SizeNumTotalNew += $SizeAll[$i][$j];
			}
		}
		my $Multi = $SizeNumTotal / $SizeNumTotalNew;
		my $tId = 0;
		open(BEFORE,"cat $File4SizeStatBeforeCorrect |") or die $! unless($File4SizeStatBeforeCorrect =~ /\.gz$/);
		open(BEFORE,"zcat $File4SizeStatBeforeCorrect |") or die $! if($File4SizeStatBeforeCorrect =~ /\.gz$/);
		open(AFTER,"> $File4CorrectedInfo") or die $! unless($File4CorrectedInfo =~ /\.gz$/);
		open(AFTER,"| gzip > $File4CorrectedInfo") or die $! if($File4CorrectedInfo =~ /\.gz$/);
		while(my $Line = <BEFORE>)
		{
			if($Line =~ /^#/)
			{
				print AFTER $Line;
				next;
			}
			chomp $Line;
			my @Cols = split /\t/, $Line;
			
			for my $i (3 .. $#Cols)
			{
				$Cols[$i] = sprintf("%.8f",$Cols[$i] - $SizeAll[$tId][$i - 3] * $Multi);
			}
			print AFTER join("\t",@Cols),"\n";
			$tId ++;
		}
		close BEFORE;
		close AFTER;
	}
}
printf "[ %s ] The end.\n",TimeString(time,$BeginTime);


######### Sub functions ##########
