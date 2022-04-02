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
my ($File4Ori,$File4Stat,$Mode,$MinDepth,$Dir);
my $HelpInfo = <<USAGE;

 $ThisScriptName
 Auther: zhangdong_xie\@foxmail.com

  This script was used to calculate the accumulated percent or pure percent from number count.
  本脚本用于计算频率或者累积频率；

 -i      ( Required ) File generated by InsertSizeInfoCollect;
 -o      ( Required ) File for the stat;

 -mode   ( Optional ) Mode (percent or accum, default: accum);
 -min    ( Optional ) Minimal Depth (default: 1);
 -bin    ( Optional ) List for searching of related bin or scripts; 
 -h      ( Optional ) Help infomation;

USAGE

GetOptions(
	'i=s' => \$File4Ori,
	'o=s' => \$File4Stat,
	'mode:s' => \$Mode,
	'min:i' => \$MinDepth,
	'bin:s' => \$BinList,
	'h!' => \$HelpFlag
) or die $HelpInfo;

if($HelpFlag || !$File4Ori || !$File4Stat)
{
	die $HelpInfo;
}
else
{
	$BeginTime = ScriptBegin(0,$ThisScriptName);
	IfFileExist($File4Ori);
	$Dir = dirname $File4Stat;
	IfDirExist($Dir);
	$MinDepth = 1 unless($MinDepth);
	
	#$BinList = BinListGet() if(!$BinList);
	#$Ref = BinSearch("Reference",$BinList);
}

if(1)
{
	open(ORI,"cat $File4Ori |") or die $! unless($File4Ori =~ /\.gz$/);
	open(ORI,"zcat $File4Ori |") or die $! if($File4Ori =~ /\.gz$/);
	open(STAT,"> $File4Stat") or die $! unless($File4Stat =~ /\.gz$/);
	open(STAT,"| gzip > $File4Stat") or die $! if($File4Stat =~ /\.gz$/);
	my $InitialId = 0;
	while(my $Line = <ORI>)
	{
		chomp $Line;
		my @Cols = split /\t/, $Line;
		if($Line =~ /^#/)
		{
			for my $i (0 .. $#Cols)
			{
				unless($Cols[$i] =~ /\D/)
				{
					$InitialId = $i;
					last;
				}
			}
			print STAT $Line,"\n";
			next;
		}
		
		my @Per = @Cols[$InitialId .. $#Cols];
		my $Total = 0;
		for my $i (0 .. $#Per)
		{
			$Total += $Per[$i];
		}
		next if($Total < $MinDepth);
		
		@Per = @{Size2Percent(\@Per)} if($Mode eq "percent");
		@Per = @{Size2AccumPercent(\@Per)} if($Mode eq "accum");
		for my $i (0 .. $#Per)
		{
			$Per[$i] = sprintf("%.8f",$Per[$i]);
		}
		print STAT join("\t",@Cols[0 .. $InitialId - 1],@Per),"\n";
	}
	close ORI;
	close STAT;
}
printf "[ %s ] The end.\n",TimeString(time,$BeginTime);


######### Sub functions ##########
