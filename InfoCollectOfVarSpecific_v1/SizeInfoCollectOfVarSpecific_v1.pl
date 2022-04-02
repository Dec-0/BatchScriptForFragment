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
my ($Bam,$VarList,$Dir,$Prefix,$Perl,$Samtools,$VarRelatedReadsExtractScript,$BamFltWithIdScript);
my $HelpInfo = <<USAGE;

 $ThisScriptName
 Auther: zhangdong_xie\@foxmail.com

  This script can be used to collect insert size info for specific variants.

  本脚本用于从bam中获取变异列表相关的reads insert size信息。

 -i      ( Required ) Bam file;
 -v      ( Requried ) Variants\' list with at least 5 coloumes 'chr from to ref A-Base ...', it could be snv or indel;
                      变异信息列表，最少5列（to坐标可以为“-”），假如同一位置有多种碱基变化，则继续在后面添加B-Base、C-Base ...
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
	IfFileExist($Bam,$VarList);
	$Dir = IfDirExist($Dir);
	
	$BinList = BinListGet() if(!$BinList);
	$Perl = BinSearch("Perl",$BinList);
	$Samtools = BinSearch("Samtools",$BinList);
	$VarRelatedReadsExtractScript = BinSearch("VarRelatedReadsExtractScript",$BinList);
	$BamFltWithIdScript = BinSearch("BamFltWithIdScript",$BinList);
}

if(1)
{
	my $BaseNum = 0;
	$BaseNum = `cat $VarList | head -n1 | sed 's/\\t/\\n/g' | wc -l` unless($VarList =~ /\.gz$/);
	$BaseNum = `zcat $VarList | head -n1 | sed 's/\\t/\\n/g' | wc -l` if($VarList =~ /\.gz$/);
	chomp $BaseNum;
	# 默认第1个碱基是第4位;
	$BaseNum -= 4;
	die "[ Error ] No column for alt base\n" unless($BaseNum >= 1);
	print "[ Info ] The number of variants' string is $BaseNum\n";
	
	# 逐个碱基列进行处理;
	for my $i (1 .. $BaseNum)
	{
		my $ColId = $i + 4;
		my $tVarList = $Dir . "/" . $Prefix . ".VarList.Var" . $i . ".xls";
		open(ORI,"cat $VarList | grep -v ^# | cut -f 1-4,$ColId |") or die $! unless($VarList =~ /\.gz$/);
		open(ORI,"zcat $VarList | grep -v ^# | cut -f 1-4,$ColId |") or die $! if($VarList =~ /\.gz$/);
		open(SPLIT,"> $tVarList") or die $!;
		while(my $Line = <ORI>)
		{
			chomp $Line;
			my @Cols = split /\t/, $Line;
			# 针对To坐标缺失的情况;
			if($Cols[3])
			{
				$Cols[2] = $Cols[1] + length($Cols[3]) - 1;
			}
			else
			{
				$Cols[2] = $Cols[1];
			}
			print SPLIT join("\t",@Cols),"\n";
		}
		close ORI;
		close SPLIT;
		
		# 搜集所有相关的reads信息，减少bam大小以提高速度;
		if(1)
		{
			`$Perl $VarRelatedReadsExtractScript -v $tVarList -b $Bam -o $Dir/$Prefix\.Var$i\.tmp.bam`;
			my %DupInfo = ();
			open(LIST,"$Samtools view $Dir/$Prefix\.Var$i\.tmp.bam | cut -f 1 |") or die $!;
			open(MLIST,"| gzip > $Dir/$Prefix\.Var$i\.ReadList.gz") or die $!;
			while(my $Line = <LIST>)
			{
				next if($DupInfo{$Line});
				print MLIST $Line;
				
				$DupInfo{$Line} = 1;
			}
			close LIST;
			close MLIST;
			`perl $BamFltWithIdScript -i $Bam -r $Dir/$Prefix\.Var$i\.ReadList.gz -o $Dir/$Prefix\.Var$i\.bam`;
			`rm $Dir/$Prefix\.Var$i\.tmp.bam`;
		}
		
		# 逐个变异统计;
		if(1)
		{
			open(SD,"> $Dir/$Prefix\.InsertSizeDistr.Var$i\.xls") or die $!;
			print SD join("\t","#Chr","From","To","Ref","Alt","Size"),"\n";
			# 为了加快分析速度，bam仍需继续拆分，一次500个变异;
			my $SPlitNum = 500;
			my $TotalVarNum = `cat $tVarList | grep -v ^# | wc -l`;
			chomp $TotalVarNum;
			my $MaxSplitNum = int($TotalVarNum / $SPlitNum);
			$MaxSplitNum -- unless($TotalVarNum % $SPlitNum > 0);
			for my $n (0 .. $MaxSplitNum)
			{
				my $MinCol = $n * $SPlitNum;
				my $MaxCol = $MinCol + $SPlitNum;
				`cat $tVarList | grep -v ^# | awk '{if(NR > $MinCol && NR <= $MaxCol){print \$0}}' > $Dir/$Prefix\.VarList.Var$i\.Split$n\.xls`;
				`$Perl $VarRelatedReadsExtractScript -v $Dir/$Prefix\.VarList.Var$i\.Split$n\.xls -b $Dir/$Prefix\.Var$i\.bam -o $Dir/$Prefix\.Var$i\.Split$n\.bam`;
				
				open(LIST,"< $Dir/$Prefix\.VarList.Var$i\.Split$n\.xls") or die $!;
				while(my $Line = <LIST>)
				{
					chomp $Line;
					`echo -e '$Line' > $Dir/$Prefix\.InsertSizeDistr.Var$i\.tmp.xls`;
					`$Perl $VarRelatedReadsExtractScript -v $Dir/$Prefix\.InsertSizeDistr.Var$i\.tmp.xls -b $Dir/$Prefix\.Var$i\.Split$n\.bam -o $Dir/$Prefix\.InsertSizeDistr.Var$i\.tmp.bam`;
					
					my %DupFlag = ();
					my @SizeList = ();
					open(SM,"$Samtools view $Dir/$Prefix\.InsertSizeDistr.Var$i\.tmp.bam | cut -f 1,9 |") or die $!;
					while(my $Line = <SM>)
					{
						chomp $Line;
						my ($Id,$Size) = split /\t/, $Line;
						next if($DupFlag{$Id});
						$DupFlag{$Id} = 1;
						
						next if($Size == 0);
						push @SizeList, abs($Size);
					}
					close SM;
					
					@SizeList = @{PureNumSort(\@SizeList)};
					my $SizeListString = join(",",@SizeList);
					$SizeListString = "-" unless($SizeListString);
					`rm $Dir/$Prefix\.InsertSizeDistr.Var$i\.tmp.xls $Dir/$Prefix\.InsertSizeDistr.Var$i\.tmp.bam`;
					print SD join("\t",$Line,$SizeListString),"\n";
				}
				close LIST;
				
				`rm $Dir/$Prefix\.VarList.Var$i\.Split$n\.xls $Dir/$Prefix\.Var$i\.Split$n\.bam`;
			}
			close SD;
		}
		`rm $tVarList`;
	}
}
printf "[ %s ] The end.\n",TimeString(time,$BeginTime);


######### Sub functions ##########
