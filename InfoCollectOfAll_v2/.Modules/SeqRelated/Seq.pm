# Package Name
package SeqRelated::Seq;

# Exported name
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(RefGet GCNum);

use FindBin qw($Bin);

sub RefGet
{
	my ($RefGen,$Chr,$From,$To,$BedtoolsBin,$OriFlag) = @_;
	
	die "[ Error ] End position can not be smaller than start position when getting seq ($Chr,$From,$To) from ref($RefGen).\n" if($From > $To);
	die "[ Error ] Reference not exist ($RefGen).\n" unless(-e $RefGen);
	$BedtoolsBin = "bedtools" unless($BedtoolsBin);
	
	# 这里默认输入的是1-based的坐标，而bedtools针对的是0-based的坐标，因此假如想取得坐标[a,b]的碱基，则需要转换成(a - 1, b);
	$From --;
	my $Ref = `echo -e '$Chr\\t$From\\t$To' | $BedtoolsBin getfasta -fi $RefGen -bed - | tail -n 1`;
	chomp $Ref;
	# 假如不需要小写变大写;
	$Ref = uc($Ref) unless($OriFlag);
	
	return $Ref;
}

sub GCNum
{
	my $String = $_[0];
	my $GCNum = 0;
	
	my @Base = split //, $String;
	for my $i (0 .. $#Base)
	{
		$GCNum ++ if($Base[$i] eq "G" || $Base[$i] eq "g" || $Base[$i] eq "C" || $Base[$i] eq "c");
	}
	
	return $GCNum;
}

1;
