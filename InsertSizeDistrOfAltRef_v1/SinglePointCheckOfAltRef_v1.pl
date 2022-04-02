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
my ($VarInfo,$File4Ref,$File4Alt,$Dir,$Prefix,$RScript,$SizeDistrVisualScript,$Script4Density);
my $HelpInfo = <<USAGE;

 $ThisScriptName
 Auther: zhangdong_xie\@foxmail.com

  This script can be used to debug for one single point.
  
  针对单个点，可视化展示Ref、Alt的片段长度密度、累积密度差异，以及EndPos的分布差异。

 -var    ( Required ) VarInfo like 'chr1,123456,A,G';
 -ref    ( Required ) InsertInfo for reference;
                      Ref碱基对应的特定插入片段信息。文件至少有6列：染色体号、起始坐标、结束坐标、ref序列、当前序列（ref或alt）、片段长度信息（逗号分隔）。
 -alt    ( Required ) InsertInfo for alt;
                      Alt碱基对应的特定插入片段信息；
 -o      ( Required ) Directory for result;
 -prefix ( Required ) Prefix of the output files;

 -bin    ( Optional ) List for searching of related bin or scripts; 
 -h      ( Optional ) Help infomation;

USAGE

GetOptions(
	'var=s' => \$VarInfo,
	'ref=s' => \$File4Ref,
	'alt=s' => \$File4Alt,
	'o=s' => \$Dir,
	'prefix=s' => \$Prefix,
	'bin:s' => \$BinList,
	'h!' => \$HelpFlag
) or die $HelpInfo;

if($HelpFlag || !$VarInfo || !$File4Ref || !$File4Alt || !$Dir || !$Prefix)
{
	die $HelpInfo;
}
else
{
	$BeginTime = ScriptBegin(0,$ThisScriptName);
	$Dir = IfDirExist($Dir);
	
	$BinList = BinListGet() if(!$BinList);
	$RScript = BinSearch("RScript",$BinList);
	$SizeDistrVisualScript = BinSearch("SizeDistrVisualScript",$BinList);
	$Script4Density = BinSearch("Script4Density",$BinList);
}

if(1)
{
	my ($Chr,$Pos,$Ref,$Alt) = split /,/, $VarInfo;
	
	my ($MinSize,$MaxSize) = (100,300);
	# Ref部分；
	my ($File4RefNum,$File4RefPer,$File4RefPerAccum,$File4RefEndPos) = ();
	if(1)
	{
		my $tRefInfo = "";
		$tRefInfo = `cat $File4Ref | awk '{if(\$1 == "$Chr" && \$2 == "$Pos" && \$4 == "$Ref" && \$5 == "$Ref"){print \$0}}' | cut -f 6,7` unless($File4Ref =~ /\.gz$/);
		$tRefInfo = `zcat $File4Ref | awk '{if(\$1 == "$Chr" && \$2 == "$Pos" && \$4 == "$Ref" && \$5 == "$Ref"){print \$0}}' | cut -f 6,7` if($File4Ref =~ /\.gz$/);
		chomp $tRefInfo;
		my ($Info4RefSize,$Info4RefPos) = split /\t/, $tRefInfo;
		$Info4RefSize = "" if($Info4RefSize eq "-");
		$Info4RefPos = "" if($Info4RefPos eq "-");
		my @RefSize = split /,/, $Info4RefSize;
		my @RefPos = split /,/, $Info4RefPos;
		
		# 片段长度部分;
		my %SizeInfo4Ref = ();
		my $TotalRef = 0;
		for my $i (0 .. $#RefSize)
		{
			next unless($RefSize[$i] >= $MinSize && $RefSize[$i] <= $MaxSize);
			$SizeInfo4Ref{$RefSize[$i]} = 0 unless($SizeInfo4Ref{$RefSize[$i]});
			$SizeInfo4Ref{$RefSize[$i]} ++;
			$TotalRef ++;
		}
		$File4RefNum = $Dir . "/" . $Prefix . ".Ref.Number.xls";
		$File4RefPer = $Dir . "/" . $Prefix . ".Ref.Percent.xls";
		$File4RefPerAccum = $Dir . "/" . $Prefix . ".Ref.PercentAccum.xls";
		open(RN,"> $File4RefNum") or die $!;
		open(RP,"> $File4RefPer") or die $!;
		open(RPA,"> $File4RefPerAccum") or die $!;
		my $tRefAccumPer = 0;
		foreach my $Size (sort {$a <=> $b} keys %SizeInfo4Ref)
		{
			print RN join("\t",$Size,$SizeInfo4Ref{$Size}),"\n";
			my $tPer = sprintf("%.5f",$SizeInfo4Ref{$Size} / $TotalRef);
			print RP join("\t",$Size,$tPer),"\n";
			
			$tRefAccumPer += $SizeInfo4Ref{$Size} / $TotalRef;
			$tPer = sprintf("%.5f",$tRefAccumPer);
			print RPA join("\t",$Size,$tPer),"\n";
		}
		close RN;
		close RP;
		close RPA;
		
		# EndPosition部分；
		$File4RefEndPos = $Dir . "/" . $Prefix . ".Ref.EndPos.xls";
		open(AP,"> $File4RefEndPos") or die $!;
		for my $i (0 .. $#RefPos)
		{
			next unless($RefSize[$i] >= $MinSize && $RefSize[$i] <= $MaxSize);
			my $LeftPos = 1 - $RefPos[$i];
			my $RightPos = $RefSize[$i] - $RefPos[$i];
			print AP join("\t","Ref",$LeftPos),"\n";
			print AP join("\t","Ref",$RightPos),"\n";
		}
		close AP;
	}
	# Alt部分；
	my ($File4AltNum,$File4AltPer,$File4AltPerAccum,$File4AltEndPos) = ();
	if(1)
	{
		my $tAltInfo = "";
		$tAltInfo = `cat $File4Alt | awk '{if(\$1 == "$Chr" && \$2 == "$Pos" && \$4 == "$Ref" && \$5 == "$Alt"){print \$0}}' | cut -f 6,7` unless($File4Alt =~ /\.gz$/);
		$tAltInfo = `zcat $File4Alt | awk '{if(\$1 == "$Chr" && \$2 == "$Pos" && \$4 == "$Ref" && \$5 == "$Alt"){print \$0}}' | cut -f 6,7` if($File4Alt =~ /\.gz$/);
		chomp $tAltInfo;
		my ($Info4AltSize,$Info4AltPos) = split /\t/, $tAltInfo;
		$Info4AltSize = "" if($Info4AltSize eq "-");
		$Info4AltPos = "" if($Info4AltPos eq "-");
		my @AltSize = split /,/, $Info4AltSize;
		my @AltPos = split /,/, $Info4AltPos;
		
		# 片段长度部分;
		my %SizeInfo4Alt = ();
		my $TotalAlt = 0;
		for my $i (0 .. $#AltSize)
		{
			next unless($AltSize[$i] >= $MinSize && $AltSize[$i] <= $MaxSize);
			$SizeInfo4Alt{$AltSize[$i]} = 0 unless($SizeInfo4Alt{$AltSize[$i]});
			$SizeInfo4Alt{$AltSize[$i]} ++;
			$TotalAlt ++;
		}
		$File4AltNum = $Dir . "/" . $Prefix . ".Alt.Number.xls";
		$File4AltPer = $Dir . "/" . $Prefix . ".Alt.Percent.xls";
		$File4AltPerAccum = $Dir . "/" . $Prefix . ".Alt.PercentAccum.xls";
		open(RN,"> $File4AltNum") or die $!;
		open(RP,"> $File4AltPer") or die $!;
		open(RPA,"> $File4AltPerAccum") or die $!;
		my $tAltAccumPer = 0;
		foreach my $Size (sort {$a <=> $b} keys %SizeInfo4Alt)
		{
			print RN join("\t",$Size,$SizeInfo4Alt{$Size}),"\n";
			my $tPer = sprintf("%.5f",$SizeInfo4Alt{$Size} / $TotalAlt);
			print RP join("\t",$Size,$tPer),"\n";
			
			$tAltAccumPer += $SizeInfo4Alt{$Size} / $TotalAlt;
			$tPer = sprintf("%.5f",$tAltAccumPer);
			print RPA join("\t",$Size,$tPer),"\n";
		}
		close RN;
		close RP;
		close RPA;
		
		# EndPosition部分；
		$File4AltEndPos = $Dir . "/" . $Prefix . ".Alt.EndPos.xls";
		open(AP,"> $File4AltEndPos") or die $!;
		for my $i (0 .. $#AltPos)
		{
			next unless($AltSize[$i] >= $MinSize && $AltSize[$i] <= $MaxSize);
			my $LeftPos = 1 - $AltPos[$i];
			my $RightPos = $AltSize[$i] - $AltPos[$i];
			print AP join("\t","Alt",$LeftPos),"\n";
			print AP join("\t","Alt",$RightPos),"\n";
		}
		close AP;
	}
	
	# 绘图；
	if(1)
	{
		# 密度；
		my $Pdf4Percent = $Dir . "/" . $Prefix . ".Percent.pdf";
		`$RScript $SizeDistrVisualScript $File4RefPer $File4AltPer $Pdf4Percent '$Chr,$Pos,$Ref,$Alt' '' '片段长度' '密度' 'Ref Base' 'Alt Base'`;
		# 累积密度；
		my $Pdf4Accum = $Dir . "/" . $Prefix . ".AccumPercent.pdf";
		`$RScript $SizeDistrVisualScript $File4RefPerAccum $File4AltPerAccum $Pdf4Accum '$Chr,$Pos,$Ref,$Alt' '' '片段长度' '累积密度' 'Ref Base' 'Alt Base'`;
		# 末端坐标密度图；
		my $Pdf4EndPos = $Dir . "/" . $Prefix . ".EndPos.pdf";
		`cat $File4RefEndPos > $Dir/$Prefix\.EndPos.xls`;
		`cat $File4AltEndPos >> $Dir/$Prefix\.EndPos.xls`;
		`$RScript $Script4Density $Pdf4EndPos $Dir/$Prefix\.EndPos.xls '$Chr,$Pos,$Ref,$Alt' '相对位置' '密度'`;
	}
}
printf "[ %s ] The end.\n",TimeString(time,$BeginTime);


######### Sub functions ##########
