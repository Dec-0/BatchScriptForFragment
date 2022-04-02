# Package Name
package MGDRelated::InsertSize;

# Exported name
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(SizeInfoGet AccumDensity AccumDensityOnHash AccumPercent AccumPercentOnHash ISNumInRange FreqInSpecificSize);
use Sort::PureNum;

# 获得对应的插入片段信息;
sub SizeInfoGet
{
	my ($File,$Chr,$Pos,$Seq) = @_;
	my @Size = ();
	
	my $InsertString = "";
	$InsertString = `cat $File | awk '{if(\$1 == "$Chr" && \$2 == "$Pos" && \$5 == "$Seq"){print \$0;exit;}}' | cut -f 6` unless($File =~ /\.gz$/);
	$InsertString = `zcat $File | awk '{if(\$1 == "$Chr" && \$2 == "$Pos" && \$5 == "$Seq"){print \$0;exit;}}' | cut -f 6` if($File =~ /\.gz$/);
	chomp $InsertString;
	die "[ Error ] No InsertInfo for $Chr,$Pos,$Alt in $File\n" unless($InsertString);
	$InsertString = "" if($InsertString eq "-");
	@Size = split /,/, $InsertString;
	
	return \@Size;
}

# 统计插入片段的累积概率密度(0~1)；
sub AccumDensity
{
	my $BMin = $_[0];
	my $BMax = $_[1];
	my $BLen = $_[2];
	my @InsertInfo = @{$_[3]};
	
	my @AccumDensity = ();
	$AccumDensity[0] = 0;
	for my $i (0 .. $#InsertInfo)
	{
		for my $j ($BMin .. $BMax)
		{
			$AccumDensity[$j - $BMin + 1] = 0 unless($AccumDensity[$j - $BMin + 1]);
			$AccumDensity[$j - $BMin + 1] ++ if($InsertInfo[$i] >= $j - $BLen && $InsertInfo[$i] <= $j + $BLen);
			last if($InsertInfo[$i] < $j - $BLen);
		}
	}
	my $AccumNum = 0;
	my $Len = $BMax - $BMin + 1;
	for my $i (0 .. $Len)
	{
		$AccumDensity[$i] = 0 unless($AccumDensity[$i]);
		$AccumNum += $AccumDensity[$i];
		$AccumDensity[$i] = $AccumNum;
	}
	if($AccumNum == 0)
	{
		@AccumDensity = ();
		return \@AccumDensity;
	}
	for my $i (0 .. $#AccumDensity)
	{
		$AccumDensity[$i] = $AccumDensity[$i] / $AccumNum;
	}
	
	return \@AccumDensity;
}
sub AccumDensityOnHash
{
	my $BMin = $_[0];
	my $BMax = $_[1];
	my $BLen = $_[2];
	my %InsertInfo = %{$_[3]};
	
	my @AccumDensity = ();
	$AccumDensity[0] = 0;
	foreach my $Size (keys %InsertInfo)
	{
		for my $j ($BMin .. $BMax)
		{
			$AccumDensity[$j - $BMin + 1] = 0 unless($AccumDensity[$j - $BMin + 1]);
			$AccumDensity[$j - $BMin + 1] += $InsertInfo{$Size} if($Size >= $j - $BLen && $Size <= $j + $BLen);
		}
	}
	my $AccumNum = 0;
	my $Len = $BMax - $BMin + 1;
	for my $i (0 .. $Len)
	{
		$AccumDensity[$i] = 0 unless($AccumDensity[$i]);
		$AccumNum += $AccumDensity[$i];
		$AccumDensity[$i] = $AccumNum;
	}
	if($AccumNum == 0)
	{
		@AccumDensity = ();
		return \@AccumDensity;
	}
	for my $i (0 .. $#AccumDensity)
	{
		$AccumDensity[$i] = $AccumDensity[$i] / $AccumNum;
	}
	
	return \@AccumDensity;
}


# 统计插入片段的累积百分比差异(0~1)；
sub AccumPercent
{
	my $BMin = $_[0];
	my $BMax = $_[1];
	my @InsertInfo = @{$_[2]};
	
	my @AccumPercent = ();
	$AccumPercent[0] = 0;
	for my $i (0 .. $#InsertInfo)
	{
		next if($InsertInfo[$i] < $BMin || $InsertInfo[$i] > $BMax);
		
		$AccumPercent[$InsertInfo[$i] - $BMin + 1] = 0 unless($AccumPercent[$InsertInfo[$i] - $BMin + 1]);
		$AccumPercent[$InsertInfo[$i] - $BMin + 1] ++;
	}
	my $AccumNum = 0;
	my $Len = $BMax - $BMin + 1;
	for my $i (0 .. $Len)
	{
		$AccumPercent[$i] = 0 unless($AccumPercent[$i]);
		$AccumNum += $AccumPercent[$i];
		$AccumPercent[$i] = $AccumNum;
	}
	if($AccumNum == 0)
	{
		@AccumPercent = ();
		return \@AccumPercent;
	}
	for my $i (0 .. $#AccumPercent)
	{
		$AccumPercent[$i] = $AccumPercent[$i] / $AccumNum;
	}
	
	return \@AccumPercent;
}
sub AccumPercentOnHash
{
	my $BMin = $_[0];
	my $BMax = $_[1];
	my %InsertInfo = %{$_[2]};
	
	my @AccumPercent = ();
	$AccumPercent[0] = 0;
	foreach my $Size (keys %InsertInfo)
	{
		next if($Size < $BMin || $Size > $BMax);
		
		$AccumPercent[$Size - $BMin + 1] = 0 unless($AccumPercent[$Size - $BMin + 1]);
		$AccumPercent[$Size - $BMin + 1] += $InsertInfo{$Size};
	}
	my $AccumNum = 0;
	my $Len = $BMax - $BMin + 1;
	for my $i (0 .. $Len)
	{
		$AccumPercent[$i] = 0 unless($AccumPercent[$i]);
		$AccumNum += $AccumPercent[$i];
		$AccumPercent[$i] = $AccumNum;
	}
	if($AccumNum == 0)
	{
		@AccumPercent = ();
		return \@AccumPercent;
	}
	for my $i (0 .. $#AccumPercent)
	{
		$AccumPercent[$i] = $AccumPercent[$i] / $AccumNum;
	}
	
	return \@AccumPercent;
}

# 用于判断特定区间内的片段数量
sub ISNumInRange
{
	my @SizeInfo = @{$_[0]};
	my $MinSize = $_[1];
	my $MaxSize = $_[2];
	
	my $Num = 0;
	for my $i (0 .. $#SizeInfo)
	{
		$Num ++ if($SizeInfo[$i] >= $MinSize && $SizeInfo[$i] <= $MaxSize);
	}
	
	return $Num;
}

sub FreqInSpecificSize
{
	my @SizeInfo4Ref = @{$_[0]};
	my @SizeInfo4Alt = @{$_[1]};
	my $MinSize = $_[2];
	my $MaxSize = $_[3];
	my $Flag4Per = $_[4];
	$Flag4Per = 0 unless($Flag4Per);
	
	if($Flag4Per)
	{
		die "[ Error ] MinSize not in range [0,1]\n" unless($MinSize >= 0 && $MinSize <= 1);
		die "[ Error ] MaxSize not in range [0,1]\n" unless($MaxSize >= 0 && $MaxSize <= 1);
		die "[ Error ] MinSize >= MaxSize\n" unless($MaxSize > $MinSize);
		my @SizeInfo4All = @SizeInfo4Ref;
		push @SizeInfo4All, @SizeInfo4Alt;
		@SizeInfo4All = @{PureNumSort(\@SizeInfo4All)};
		my $TotalNum = @SizeInfo4All;
		my $MinId = int($TotalNum * $MinSize + 0.5);
		my $MaxId = int($TotalNum * $MaxSize + 0.5);
		$MinId -- if($MinId > 0);
		#$MaxId -- if($MaxId > 0);
		$MinSize = $SizeInfo4All[$MinId];
		$MaxSize = $SizeInfo4All[$MaxId];
	}
	
	my $RefNum = &ISNumInRange(\@SizeInfo4Ref,$MinSize,$MaxSize);
	my $AltNum = &ISNumInRange(\@SizeInfo4Alt,$MinSize,$MaxSize);
	return "-" unless($RefNum + $AltNum > 0);
	my $Freq = $AltNum / ($RefNum + $AltNum);
	
	return $Freq;
}

1;
