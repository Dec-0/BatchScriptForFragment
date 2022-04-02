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
my ($Bam,$VarList,$Dir,$Prefix);
my $Samtools;
my $HelpInfo = <<USAGE;

 $ThisScriptName
 Auther: zhangdong_xie\@foxmail.com

  This script can be used to collect insert size (Related but not cover) info for specific variants.
  
  本脚本用于收集变异相关但Reads并不覆盖的相关插入片段信息。

 -i      ( Required ) Bam file;
                      Bam必须要按坐标排序后的，且需要建好index;
 -v      ( Required ) A list which records the variants\' info (one snv/indel a line, include \'chr,from,to,ref,alt\' at least);
 -o      ( Requried ) Directory for result logging;
 -prefix ( Requried ) Prefix of the logging file;

 -bin    ( Optional ) List for searching of related bin or scripts; 
 -h      ( Optional ) Help infomation;

USAGE

GetOptions(
	'i=s' => \$Bam,
	'v=s' => \$VarList,
	'o=s' => \$Dir,
	'prefix=s' => \$Prefix,
	'bin:s' => \$BinList,
	'h!' => \$HelpFlag
) or die $HelpInfo;

if($HelpFlag || !$Bam || !$VarList || !$Dir || !$Prefix)
{
	die $HelpInfo;
}
else
{
	$BeginTime = ScriptBegin(0,$ThisScriptName);
	
	die "[ Error ] Not bam file ($Bam).\n" unless($Bam =~ /\.bam$/);
	IfFileExist($Bam,$VarList);
	$Dir = IfDirExist($Dir);
	
	$BinList = BinListGet() if(!$BinList);
	$Samtools = BinSearch("Samtools",$BinList);
}

if(1)
{
	# 确定Reads读长;
	my $ReadLen = 0;
	if(1)
	{
		my $MaxId = 1000;
		my $tId = 1;
		open(TBAM,"$Samtools view -F 0xF04 $Bam | cut -f 10 |") or die $!;
		while(my $Line = <TBAM>)
		{
			my $tLen = length($Line);
			$ReadLen = $tLen if($tLen > $ReadLen);
			
			last if($tId > $MaxId);
			$tId ++;
		}
		close TBAM;
		$ReadLen --;
		die "[ Error ] Unproper read length $ReadLen in the first $MaxId reads.\n" unless($ReadLen > 10);
		print "[ Info ] Read length $ReadLen\n";
	}
	
	# 检查是否存在bai文件;
	if(1)
	{
		my $Bai = $Bam . ".bai";
		unless(-s $Bai)
		{
			print "[ Info ] Begin index bam file ($Bam).\n";
			`$Samtools index $Bam`;
			print "[ Info ] Bam index done.\n";
		}
	}
	
=pod
	# 创建一个临时bed，拓宽以覆盖相关没有Cover的reads信息;
	my $EBed = $Dir . "/" . $Prefix . ".extend." . RandString(10) . ".bed";
	my %CoordInfo = ();
	if(1)
	{
		my $ExtendSize = 500;
		my (@Chr,@Form,@To,@Other) = ();
		open(VARL,"zcat $VarList | grep -v ^# | cut -f 1-3 |") or die $! if($VarList =~ /\.gz$/);
		open(VARL,"cat $VarList | grep -v ^# | cut -f 1-3 |") or die $! unless($VarList =~ /\.gz$/);
		while(my $Line = <VARL>)
		{
			chomp $Line;
			my @Cols = split /\t/, $Line;
			push @Chr, $Cols[0];
			push @From, $Cols[1];
			push @To, $Cols[2];
			push @Other, "";
		}
		close VARL;
		my @tRef = ChrPosAndOther(\@Chr,\@From,\@To,\@Other);
		@Chr = @{$tRef[0]};
		@From = @{$tRef[1]};
		@To = @{$tRef[2]};
		
		open(TMP,"| $Bedtools merge -i - > $EBed") or die $!;
		for my $i (0 .. $#Chr)
		{
			push @{$CoordInfo{$Chr[$i]}{"From"}}, $From[$i];
			push @{$CoordInfo{$Chr[$i]}{"To"}}, $To[$i];
			my $tFrom = $From[$i] - 1;
			my $tTo = $To[$i] - 1;
			print TMP join("\t",$Chr[$i],$tFrom,$tTo),"\n";
		}
		close TMP;
	}
=cut
	
	# 统计信息;
	my $ExtendSize = 500;
	open(VARL,"zcat $VarList | grep -v ^# |") or die $! if($VarList =~ /\.gz$/);
	open(VARL,"cat $VarList | grep -v ^# |") or die $! unless($VarList =~ /\.gz$/);
	open(SD,"> $Dir/$Prefix\.InsertSizeRelatedButNotCover.xls") or die $!;
	print SD join("\t","#Chr","From","To","Size","DistanceFrom5End"),"\n";
	while(my $VarLine = <VARL>)
	{
		chomp $VarLine;
		my ($Chr,$From,$To,$Ref,$Other) = split /\t/, $VarLine;
		if($To eq ".")
		{
			$To = $From + length($Ref) - 1 if($Ref ne "-");
			$To = $From if($Ref eq "-");
		}
		my $EFrom = $From - $ExtendSize;
		my $ETo = $To + $ExtendSize;
		
		my (@SizeList,@LenList) = ();
		open(EBAM,"$Samtools view -F 0xF04 $Bam $Chr\:$EFrom\-$ETo | cut -f 1,4,9 | awk '{if(\$3 != 0){print \$0}}' |") or die $!;
		while(my $BamLine = <EBAM>)
		{
			chomp $BamLine;
			my ($ReadId,$Pos,$Len) = split /\t/, $BamLine;
			my ($SizeFrom,$SizeTo) = ($Pos,$Pos);
			$SizeTo = $Pos + $Len - 1 if($Len > 0);
			$SizeFrom = $Pos - $Len - 1 if($Len < 0);
			next unless($From > $SizeFrom && $To < $SizeTo);
			next if($From < $SizeFrom + $ReadLen);
			next if($To > $SizeTo - $ReadLen);
			
			push @SizeList, abs($Len);
			push @LenList, $From - $SizeFrom + 1;
		}
		close EBAM;
		
		my @tRef = NumStringSort(\@SizeList,\@LenList);
		@SizeList = @{$tRef[0]};
		@LenList = @{$tRef[1]};
		my $SizeListString = join(",",@SizeList);
		$SizeListString = "-" unless($SizeListString);
		my $LenListString = join(",",@LenList);
		$LenListString = "-" unless($LenListString);
		print SD join("\t",$Chr,$From,$To,$SizeListString,$LenListString),"\n";
	}
	close VARL;
}
printf "[ %s ] The end.\n",TimeString(time,$BeginTime);


######### Sub functions ##########
