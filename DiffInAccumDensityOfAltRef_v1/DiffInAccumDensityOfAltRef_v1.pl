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

my ($HelpFlag,$BinList,$BeginTime);
my $ThisScriptName = basename $0;
my ($File4VarList,$File4Ref,$File4Alt,$File4Result,$BMin,$BMax,$BLen,$Flag4Density,$Flag4Smooth,$Flag4Merge);
my ($Python3,$Script4Smooth,$Dir);
my $HelpInfo = <<USAGE;

 $ThisScriptName
 Auther: zhangdong_xie\@foxmail.com

  This script was used to calculate the difference of accumulated density of InsertSize from Alt and Ref.
  
  本脚本用于计算Alt和Ref对应插入片段累积密度分布的差值分布，并视参数进行平滑处理。
  
  一般情况下每个长度对应的数值就是它的数量，但是在密度模式下该长度对应的数值为以它为中心一定范围内数量的均值，以-density参数为准。
 
 -v       ( Required ) List for variants;
                       至少含有5列：染色体号、起始坐标、结束坐标、ref、alt。
 -ref     ( Required ) InsertInfo for reference;
                       Ref碱基对应的特定插入片段信息。文件至少有6列：染色体号、起始坐标、结束坐标、ref序列、当前序列（ref或alt）、片段长度信息（逗号分隔）。
 -alt     ( Required ) InsertInfo for alt;
                       Alt碱基对应的特定插入片段信息；
 -o       ( Required ) File for result;

 # 长度统计范围参数
 -min     ( Optional ) The minimal fragmental size (default: 120);
 -max     ( Optional ) The maximal fragmental size (default: 200);
 
 # 密度统计模式相关参数
 -density ( Optional ) If replace the number corresponding to one point with the mean value of adjacent points (with -density only);
                       假如用以该点为中心的连续多个点的平均值代替该点的数量值。
 -e       ( Optional ) The half extended length for mean calculation (default: 10);
 
 # 平滑处理参数
 -s       ( Optional ) If do the smooth step;
 
 # 多个点的合并处理
 -a       ( Optional ) If treat all the variants as a whole one;
                       在隐性单基因病检测项目中，可用于合并多个分型相同位点的数据。
 
 -bin     ( Optional ) List for searching of related bin or scripts; 
 -h       ( Optional ) Help infomation;

USAGE

GetOptions(
	'v=s' => \$File4VarList,
	'ref=s' => \$File4Ref,
	'alt=s' => \$File4Alt,
	'o=s' => \$File4Result,
	'min:i' => \$BMin,
	'max:i' => \$BMax,
	'density!' => \$Flag4Density,
	'e:i' => \$BLen,
	's!' => \$Flag4Smooth,
	'a!' => \$Flag4Merge,
	'bin:s' => \$BinList,
	'h!' => \$HelpFlag
) or die $HelpInfo;

if($HelpFlag || !$File4VarList || !$File4Ref || !$File4Alt || !$File4Result)
{
	die $HelpInfo;
}
else
{
	$BeginTime = ScriptBegin(0,$ThisScriptName);
	IfFileExist($File4VarList,$File4Ref,$File4Alt);
	$Dir = dirname $File4Result;
	$Dir = IfDirExist($Dir);
	
	$BinList = BinListGet() if(!$BinList);
	$Script4Smooth = BinSearch("Script4Smooth",$BinList) if($Flag4Smooth);
	$Python3 = BinSearch("Python3",$BinList) if($Flag4Smooth);
	
	$BMin = 120 unless($BMin);
	$BMax = 200 unless($BMax);
	$BLen = 10 unless($BLen);
}

if(1)
{
	# 收集InsertSize信息;
	my @Id = ();
	push @Id, "Pos_0";
	for my $i ($BMin .. $BMax)
	{
		$Id[$i - $BMin + 1] = "Pos_" . $i;
	}
	my $Header = "";
	$Header = `cat $File4VarList | grep ^# | tail -n1` unless($File4VarList =~ /\.gz$/);
	$Header = `zcat $File4VarList | grep ^# | tail -n1` if($File4VarList =~ /\.gz$/);
	chomp $Header;
	my $FirstLine = "";
	$FirstLine = `cat $File4VarList | grep -v ^# | head -n1` unless($File4VarList =~ /\.gz$/);
	$FirstLine = `zcat $File4VarList | grep -v ^# | head -n1` if($File4VarList =~ /\.gz$/);
	chomp $FirstLine;
	my @FirstCol = split /\t/, $FirstLine;
	
	open(VL,"cat $File4VarList | grep -v ^# |") or die $! unless($File4VarList =~ /\.gz$/);
	open(VL,"zcat $File4VarList | grep -v ^# |") or die $! if($File4VarList =~ /\.gz$/);
	open(AD,"> $File4Result") or die $!;
	if($Header)
	{
		print AD join("\t",$Header,@Id),"\n";
	}
	else
	{
		my @Id4Other = ();
		for my $i (5 .. $#FirstCol)
		{
			push @Id4Other, "-";
		}
		print AD join("\t","#Chr","From","To","Ref","Alt",@Id4Other,@Id),"\n";
	}
	
	my (%InsertInfo4Ref,%InsertInfo4Alt) = ();
	while(my $Line = <VL>)
	{
		chomp $Line;
		my ($Chr,$From,$To,$Ref,$Alt,$Other) = split /\t/, $Line;
		
		my @RefInsert = @{SizeInfoGet($File4Ref,$Chr,$From,$Ref)};
		my @AltInsert = @{SizeInfoGet($File4Alt,$Chr,$From,$Alt)};
		
		# 假如需要多点合并处理；
		if($Flag4Merge)
		{
			for my $i (0 .. $#RefInsert)
			{
				$InsertInfo4Ref{$RefInsert[$i]} = 0 unless($InsertInfo4Ref{$RefInsert[$i]});
				$InsertInfo4Ref{$RefInsert[$i]} ++;
			}
			for my $i (0 .. $#AltInsert)
			{
				$InsertInfo4Alt{$AltInsert[$i]} = 0 unless($InsertInfo4Alt{$AltInsert[$i]});
				$InsertInfo4Alt{$AltInsert[$i]} ++;
			}
			
			next;
		}
		
		# ----------------------------
		my (@AccumRef,@AccumAlt) = ();
		if($Flag4Density)
		{
			@AccumRef = @{AccumDensity($BMin,$BMax,$BLen,\@RefInsert)};
			@AccumAlt = @{AccumDensity($BMin,$BMax,$BLen,\@AltInsert)};
		}
		else
		{
			@AccumRef = @{AccumPercent($BMin,$BMax,\@RefInsert)};
			@AccumAlt = @{AccumPercent($BMin,$BMax,\@AltInsert)};
		}
		next unless(@AccumRef && @AccumAlt);
		
		my @DiffAccum = ();
		for my $i (0 .. $#AccumRef)
		{
			push @DiffAccum, sprintf("%.8f",$AccumAlt[$i] - $AccumRef[$i]);
		}
		
		print AD join("\t",$Line,@DiffAccum),"\n";
	}
	if($Flag4Merge)
	{
		my (@AccumRef,@AccumAlt) = ();
		if($Flag4Density)
		{
			@AccumRef = @{AccumDensityOnHash($BMin,$BMax,$BLen,\%InsertInfo4Ref)};
			@AccumAlt = @{AccumDensityOnHash($BMin,$BMax,$BLen,\%InsertInfo4Alt)};
		}
		else
		{
			@AccumRef = @{AccumPercentOnHash($BMin,$BMax,\%InsertInfo4Ref)};
			@AccumAlt = @{AccumPercentOnHash($BMin,$BMax,\%InsertInfo4Alt)};
		}
		
		my @DiffAccum = ();
		for my $i (0 .. $#AccumRef)
		{
			push @DiffAccum, sprintf("%.8f",$AccumAlt[$i] - $AccumRef[$i]);
		}
		my $PreNum = @FirstCol;
		my @PreItems = ("-") x $PreNum;
		print AD join("\t",@PreItems,@DiffAccum),"\n";
	}
	close VL;
	close AD;
	
	if($Flag4Smooth)
	{
		my $NumA = @FirstCol;
		$NumA ++;
		my $NumB = $NumA + $BMax - $BMin + 1;
		my $tFile = $File4Result . ".tmp";
		my $Return = `$Python3 $Script4Smooth -i $File4Result -o $tFile -a $NumA -b $NumB`;
		chomp $Return;
		print $Return,"\n" if($Return);
		`mv $tFile $File4Result`;
	}
}
printf "[ %s ] The end.\n",TimeString(time,$BeginTime);


######### Sub functions ##########