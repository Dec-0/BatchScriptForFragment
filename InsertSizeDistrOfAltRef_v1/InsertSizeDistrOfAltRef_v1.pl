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
my ($File4GType,$File4Ref,$File4Alt,$Dir,$Prefix,$AllFlag,$ColId4MGType,$ColId4FGType,$Flag4Rm,$MinSize,$MaxSize,$RScript,$SizeDistrVisualScript);
my $HelpInfo = <<USAGE;

 $ThisScriptName
 Auther: zhangdong_xie\@foxmail.com

  This script can be used to visualize the distribution of Ref and Alt based on fetal genotype.
  
  本脚本用于（按子代基因型分类）可视化展示Ref和Alt碱基对应插入片段的密度、及累积密度曲线。
 
 -g      ( Required ) File generaged by GenTypePairConfirm_v1;
                      需要至少5列：染色体号、起始坐标、结束坐标、Ref、Alt。假如需要过滤母本杂合、以及区分子代基因型，则需要单独指定母子基因型对应的列号；
 -ref    ( Required ) InsertInfo for reference;
                      Ref碱基对应的特定插入片段信息。文件至少有6列：染色体号、起始坐标、结束坐标、ref序列、当前序列（ref或alt）、片段长度信息（逗号分隔）。
 -alt    ( Required ) InsertInfo for alt;
                      Alt碱基对应的特定插入片段信息；
 -o      ( Required ) Directory for result;
 -prefix ( Required ) Prefix of the output files;

 # 基因型相关参数
 -m      ( Optional ) Column number for maternal genotype;
                      用于过滤得到母本基因型为AB的位点，假如不指定那么所有位点都会囊括进来。
 -f      ( Optional ) Column number for fetal genotype;
                      用于判断胎儿基因型，按胎儿基因型分类数据。
 -a      ( Optional ) If do not seperate by GenType;
                      假如不分基因型混合处理
 
 # 插入片段长度统计范围
 -min    ( Optional ) Minimal Size (50bp);
 -max    ( Optional ) Maximal Size (300bp);
 
 -r      ( Optional ) If remove the redundant files;
 -bin    ( Optional ) List for searching of related bin or scripts;
 -h      ( Optional ) Help infomation;

USAGE

GetOptions(
	'g=s' => \$File4GType,
	'ref=s' => \$File4Ref,
	'alt=s' => \$File4Alt,
	'o=s' => \$Dir,
	'prefix=s' => \$Prefix,
	'm:i' => \$ColId4MGType,
	'f:i' => \$ColId4FGType,
	'a!' => \$AllFlag,
	'r!' => \$Flag4Rm,
	'min:i' => \$MinSize,
	'max:i' => \$MaxSize,
	'bin:s' => \$BinList,
	'h!' => \$HelpFlag
) or die $HelpInfo;

if($HelpFlag || !$File4GType || !$File4Ref || !$File4Alt || !$Prefix)
{
	die $HelpInfo;
}
else
{
	$BeginTime = ScriptBegin(0,$ThisScriptName);
	IfFileExist($File4GType,$File4Ref,$File4Alt);
	$Dir = IfDirExist($Dir);
	
	$BinList = BinListGet() if(!$BinList);
	$MinSize = BinSearch("MinSize",$BinList,1) unless($MinSize);
	$MaxSize = BinSearch("MaxSize",$BinList,1) unless($MaxSize);
	$RScript = BinSearch("RScript",$BinList);
	$SizeDistrVisualScript = BinSearch("SizeDistrVisualScript",$BinList);
	print "[ Info ] Size limitation: $MinSize ~ $MaxSize\n";
}

if(1)
{
	# 按胎儿基因型分类收集插入片段长度相关信息;
	my (%RefInfo,%AltInfo,%RefNum,%AltNum) = ();
	open(GT,"cat $File4GType | grep -v ^# |") unless($File4GType =~ /\.gz$/);
	open(GT,"zcat $File4GType | grep -v ^# |") if($File4GType =~ /\.gz$/);
	while(my $Line = <GT>)
	{
		chomp $Line;
		my @Cols = split /\t/, $Line;
		# 假如存在母亲基因型，则只处理母亲杂合的位点；
		my $ColNum = @Cols;
		if($ColId4MGType && $ColId4MGType <= $ColNum)
		{
			next unless($Cols[$ColId4MGType - 1] eq "AB");
		}
		my $FGType = "None";
		$FGType = $Cols[$ColId4FGType - 1] if($ColId4FGType && $ColId4FGType <= $ColNum);
		
		my $Return = `cat $File4Ref | grep ^'$Cols[0]'\$'\\t'  | awk '{if(\$1 == "$Cols[0]" && \$2 == "$Cols[1]" && \$4 == "$Cols[3]" && \$5 == "$Cols[3]"){print \$6;exit;}}'`;
		chomp $Return;
		$Return = "" if($Return eq "-");
		my @Size = split /,/, $Return;
		for my $i (0 .. $#Size)
		{
			next unless($Size[$i] >= $MinSize && $Size[$i] <= $MaxSize);
			$RefInfo{$FGType}{$Size[$i]} = 0 unless($RefInfo{$FGType}{$Size[$i]});
			$RefInfo{$FGType}{$Size[$i]} ++;
			
			$RefNum{$FGType} = 0 unless($RefNum{$FGType});
			$RefNum{$FGType} ++;
		}
		
		$Return = `cat $File4Alt | grep ^'$Cols[0]'\$'\\t' | awk '{if(\$1 == "$Cols[0]" && \$2 == "$Cols[1]" && \$4 == "$Cols[3]" && \$5 == "$Cols[4]"){print \$6;exit;}}'`;
		chomp $Return;
		$Return = "" if($Return eq "-");
		@Size = split /,/, $Return;
		for my $i (0 .. $#Size)
		{
			next unless($Size[$i] >= $MinSize && $Size[$i] <= $MaxSize);
			$AltInfo{$FGType}{$Size[$i]} = 0 unless($AltInfo{$FGType}{$Size[$i]});
			$AltInfo{$FGType}{$Size[$i]} ++;
			
			$AltNum{$FGType} = 0 unless($AltNum{$FGType});
			$AltNum{$FGType} ++;
		}
	}
	close GT;
	
	if($AllFlag)
	{
		# 所有基因型混合处理;
		my %SizeInfoRef = ();
		my $TotalRef = 0;
		foreach my $GType (keys %RefInfo)
		{
			foreach my $Size (sort {$a <=> $b} keys %{$RefInfo{$GType}})
			{
				$SizeInfoRef{$Size} = 0 unless($SizeInfoRef{$Size});
				$SizeInfoRef{$Size} += $RefInfo{$GType}{$Size};
			}
			$TotalRef += $RefNum{$GType};
		}
		my $tFile4RefNum = $Dir . "/" . $Prefix . ".Ref.Number.xls";
		my $tFile4RefPer = $Dir . "/" . $Prefix . ".Ref.Percent.xls";
		my $tFile4RefPerAccum = $Dir . "/" . $Prefix . ".Ref.PercentAccum.xls";
		open(RN,"> $tFile4RefNum") or die $!;
		open(RP,"> $tFile4RefPer") or die $!;
		open(RPA,"> $tFile4RefPerAccum") or die $!;
		my $tRefAccumPer = 0;
		foreach my $Size (sort {$a <=> $b} keys %SizeInfoRef)
		{
			print RN join("\t",$Size,$SizeInfoRef{$Size}),"\n";
			my $tPer = sprintf("%.5f",$SizeInfoRef{$Size} / $TotalRef);
			print RP join("\t",$Size,$tPer),"\n";
			
			$tRefAccumPer += $SizeInfoRef{$Size} / $TotalRef;
			$tPer = sprintf("%.5f",$tRefAccumPer);
			print RPA join("\t",$Size,$tPer),"\n";
		}
		close RN;
		close RP;
		close RPA;
		
		
		my %SizeInfoAlt = ();
		my $TotalAlt = 0;
		foreach my $GType (keys %AltInfo)
		{
			foreach my $Size (sort {$a <=> $b} keys %{$AltInfo{$GType}})
			{
				$SizeInfoAlt{$Size} = 0 unless($SizeInfoAlt{$Size});
				$SizeInfoAlt{$Size} += $AltInfo{$GType}{$Size};
			}
			$TotalAlt += $AltNum{$GType};
		}
		my $tFile4AltNum = $Dir . "/" . $Prefix . ".Alt.Number.xls";
		my $tFile4AltPer = $Dir . "/" . $Prefix . ".Alt.Percent.xls";
		my $tFile4AltPerAccum = $Dir . "/" . $Prefix . ".Alt.PercentAccum.xls";
		open(RN,"> $tFile4AltNum") or die $!;
		open(RP,"> $tFile4AltPer") or die $!;
		open(RPA,"> $tFile4AltPerAccum") or die $!;
		my $tAltAccumPer = 0;
		foreach my $Size (sort {$a <=> $b} keys %SizeInfoAlt)
		{
			print RN join("\t",$Size,$SizeInfoAlt{$Size}),"\n";
			my $tPer = sprintf("%.5f",$SizeInfoAlt{$Size} / $TotalAlt);
			print RP join("\t",$Size,$tPer),"\n";
			
			$tAltAccumPer += $SizeInfoAlt{$Size} / $TotalAlt;
			$tPer = sprintf("%.5f",$tAltAccumPer);
			print RPA join("\t",$Size,$tPer),"\n";
		}
		close RN;
		close RP;
		close RPA;
		
		my $tPdf = $Dir . "/" . $Prefix . ".pdf";
		`$RScript $SizeDistrVisualScript $tFile4RefPer $tFile4AltPer $tPdf '子代基因型: All' '' '' '' 'Ref Base' 'Alt Base'`;
		$tPdf = $Dir . "/" . $Prefix . ".Accum.pdf";
		`$RScript $SizeDistrVisualScript $tFile4RefPerAccum $tFile4AltPerAccum $tPdf '子代基因型: All' '' '' '' 'Ref Base' 'Alt Base'`;
		
		if($Flag4Rm)
		{
			`rm $tFile4RefNum` if(-s $tFile4RefNum);
			`rm $tFile4RefPer` if(-s $tFile4RefPer);
			`rm $tFile4RefPerAccum` if(-s $tFile4RefPerAccum);
			
			`rm $tFile4AltNum` if(-s $tFile4AltNum);
			`rm $tFile4AltPer` if(-s $tFile4AltPer);
			`rm $tFile4AltPerAccum` if(-s $tFile4AltPerAccum);
		}
	}
	else
	{
		# 基因型逐个处理；
		foreach my $GType (keys %RefInfo)
		{
			print "[ Info ] Processing $GType\n";
			
			my $tRefAccumPer = 0;
			my $tFile4RefNum = $Dir . "/" . $Prefix . "." . $GType . ".Ref.Number.xls";
			my $tFile4RefPer = $Dir . "/" . $Prefix . "." . $GType . ".Ref.Percent.xls";
			my $tFile4RefPerAccum = $Dir . "/" . $Prefix . "." . $GType . ".Ref.PercentAccum.xls";
			open(RN,"> $tFile4RefNum") or die $!;
			open(RP,"> $tFile4RefPer") or die $!;
			open(RPA,"> $tFile4RefPerAccum") or die $!;
			foreach my $Size (sort {$a <=> $b} keys %{$RefInfo{$GType}})
			{
				print RN join("\t",$Size,$RefInfo{$GType}{$Size}),"\n";
				my $tPer = sprintf("%.5f",$RefInfo{$GType}{$Size} / $RefNum{$GType});
				print RP join("\t",$Size,$tPer),"\n";
				
				$tRefAccumPer += $RefInfo{$GType}{$Size} / $RefNum{$GType};
				$tPer = sprintf("%.5f",$tRefAccumPer);
				print RPA join("\t",$Size,$tPer),"\n";
			}
			close RN;
			close RP;
			close RPA;
			
			my $tAltAccumPer = 0;
			my $tFile4AltNum = $Dir . "/" . $Prefix . "." . $GType . ".Alt.Number.xls";
			my $tFile4AltPer = $Dir . "/" . $Prefix . "." . $GType . ".Alt.Percent.xls";
			my $tFile4AltPerAccum = $Dir . "/" . $Prefix . "." . $GType . ".Alt.PercentAccum.xls";
			open(RN,"> $tFile4AltNum") or die $!;
			open(RP,"> $tFile4AltPer") or die $!;
			open(RPA,"> $tFile4AltPerAccum") or die $!;
			foreach my $Size (sort {$a <=> $b} keys %{$AltInfo{$GType}})
			{
				print RN join("\t",$Size,$AltInfo{$GType}{$Size}),"\n";
				my $tPer = sprintf("%.5f",$AltInfo{$GType}{$Size} / $AltNum{$GType});
				print RP join("\t",$Size,$tPer),"\n";
				
				$tAltAccumPer += $AltInfo{$GType}{$Size} / $AltNum{$GType};
				$tPer = sprintf("%.5f",$tAltAccumPer);
				print RPA join("\t",$Size,$tPer),"\n";
			}
			close RN;
			close RP;
			close RPA;
			
			my $tPdf = $Dir . "/" . $Prefix . "." . $GType . ".pdf";
			`$RScript $SizeDistrVisualScript $tFile4RefPer $tFile4AltPer $tPdf '子代基因型: $GType' '' '' '' 'Ref Base' 'Alt Base'`;
			$tPdf = $Dir . "/" . $Prefix . "." . $GType . ".Accum.pdf";
			`$RScript $SizeDistrVisualScript $tFile4RefPerAccum $tFile4AltPerAccum $tPdf '子代基因型: $GType' '' '' '' 'Ref Base' 'Alt Base'`;
			
			if($Flag4Rm)
			{
				`rm $tFile4RefNum` if(-s $tFile4RefNum);
				`rm $tFile4RefPer` if(-s $tFile4RefPer);
				`rm $tFile4RefPerAccum` if(-s $tFile4RefPerAccum);
				
				`rm $tFile4AltNum` if(-s $tFile4AltNum);
				`rm $tFile4AltPer` if(-s $tFile4AltPer);
				`rm $tFile4AltPerAccum` if(-s $tFile4AltPerAccum);
			}
		}
	}
}
printf "[ %s ] The end.\n",TimeString(time,$BeginTime);


######### Sub functions ##########