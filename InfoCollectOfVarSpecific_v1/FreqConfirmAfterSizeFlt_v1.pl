#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
Getopt::Long::Configure qw(no_ignore_case);
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin/.Modules";
use Parameter::BinList;
use Sort::PureNum;

my ($HelpFlag,$BinList,$BeginTime);
my $ThisScriptName = basename $0;
my ($ISFile,$LogFile,$Min,$Max,$Dir,$Bedtools);
my $HelpInfo = <<USAGE;

 $ThisScriptName
 Auther: zhangdong_xie\@foxmail.com

  This script can be used to visualize the freq distri from the file generate from VarSpecInsertSizeInfo.
  
  本脚本用于根据过滤后的片段分布和数量，来重新确定变异频率。

 -i      ( Required ) File recording the insert size of A-Base and B-Base;
 -o      ( Requried ) Loggging file;

 # 片段长度范围；
 -min    ( Optional ) The minimal level for insert size;
                      整数时直接当阈值，小数时当百分比。
 -max    ( Optional ) The maximal level for insert size;
 
 -bin    ( Optional ) List for searching of related bin or scripts; 
 -h      ( Optional ) Help infomation;

USAGE

GetOptions(
	'i=s' => \$ISFile,
	'o=s' => \$LogFile,
	'min:s' => \$Min,
	'max:s' => \$Max,
	'bin:s' => \$BinList,
	'h!' => \$HelpFlag
) or die $HelpInfo;

if($HelpFlag || !$ISFile || !$LogFile)
{
	die $HelpInfo;
}
else
{
	$BeginTime = ScriptBegin(0,$ThisScriptName);
	IfFileExist($ISFile);
	$Dir = dirname $LogFile;
	$Dir = IfDirExist($Dir);
	if(defined $Min && defined $Max)
	{
		die "[ Error ] Min not pure number.\n" if($Min =~ /[^\d\.]/);
		die "[ Error ] Max not pure number.\n" if($Max =~ /[^\d\.]/);
		die "[ Error ] Min > Max.\n" if($Min > $Max);
		
		if($Min <= 1 && $Max <= 1)
		{
			die "[ Error ] Min not in [0,1]\n" if($Min < 0 || $Min > 1);
			die "[ Error ] Max not in [0,1]\n" if($Max < 0 || $Max > 1);
		}
		else
		{
			$Min = int($Min);
			$Max = int($Max);
		}
		print "[ Info ] Level: $Min ~ $Max\n";
	}
	
	$BinList = BinListGet() if(!$BinList);
	$Bedtools = BinSearch("Bedtools",$BinList);
}

if(1)
{
	# 排序;
	my $tFile = $LogFile . ".tmpSort.xls";
	if(1)
	{
		`cat $ISFile | grep ^# > $tFile` unless($ISFile =~ /\.gz$/);
		`zcat $ISFile | grep ^# > $tFile` if($ISFile =~ /\.gz$/);
		
		`cat $ISFile | grep -v ^# | $Bedtools sort -i - >> $tFile` unless($ISFile =~ /\.gz$/);
		`zcat $ISFile | grep -v ^# | $Bedtools sort -i - >> $tFile` if($ISFile =~ /\.gz$/);
	}
	
	# 转换;
	if(1)
	{
		open(SORT,"cat $tFile | grep -v ^# |") or die $!;
		open(TRANS,"> $LogFile") or die $!;
		print TRANS join("\t","#Chr","From","To","Ref","Alt","Freq","FullDepth","AltDepth"),"\n";
		while(my $Line = <SORT>)
		{
			chomp $Line;
			my ($Chr,$From,$To,$Ref,$Alt,$RefString,$AltString) = split /\t/, $Line;
			$RefString = "" if($RefString eq "-");
			$AltString = "" if($AltString eq "-");
			
			# 过滤;
			my (@Ref,@Alt) = ();
			@Ref = split /,/, $RefString;
			@Alt = split /,/, $AltString;
			# 排除空值;
			next unless(@Ref || @Alt);
			if(defined $Min && defined $Max)
			{
				my ($FinalMin,$FinalMax) = ($Min,$Max);
				if($Min <= 1 && $Max <= 1)
				{
					my @Size = @Ref;
					push @Size, @Alt;
					@Size = @{PureNumSort(\@Size)};
					my $SizeNum = @Size;
					$FinalMin = $Size[int($SizeNum * $Min + 0.5)];
					$FinalMax = $Size[int($SizeNum * $Max + 0.5)];
					print "[ Info ] Transed Level: $FinalMin ~ $FinalMax\n";
				}
				
				my (@tRef,@tAlt) = ();
				for my $i (0 .. $#Ref)
				{
					next if($Ref[$i] < $FinalMin || $Ref[$i] > $FinalMax);
					push @tRef, $Ref[$i];
				}
				for my $i (0 .. $#Alt)
				{
					next if($Alt[$i] < $FinalMin || $Alt[$i] > $FinalMax);
					push @tAlt, $Alt[$i];
				}
				@Ref = @tRef;
				@Alt = @tAlt;
			}
			
			# 统计频率;
			my $RefDp = @Ref;
			my $AltDp = @Alt;
			my $FullDp = $RefDp + $AltDp;
			# 排除空值;
			next unless($FullDp);
			my $Freq = sprintf("%.4f",$AltDp / $FullDp);
			print TRANS join("\t",$Chr,$From,$To,$Ref,$Alt,$Freq,$FullDp,$AltDp),"\n";
		}
		close SORT;
		close TRANS;
	}
	`rm $tFile`;
}
printf "[ %s ] The end.\n",TimeString(time,$BeginTime);


######### Sub functions ##########
