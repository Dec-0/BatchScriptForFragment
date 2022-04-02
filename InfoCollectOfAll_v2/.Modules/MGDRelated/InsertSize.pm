# Package Name
package MGDRelated::InsertSize;

# Exported name
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(SizeInfoGet SizeInfoGetOnBed SizeInfoGetOnArea SizeInfoGetOnChr SizeInfoOfSingleOnBed SizeInfoOfMultiArea SizeCountMerge SizeCountBasedOnLo2014 SubtractOfDistrMax SubtractOfDistrMean Size2Percent Size2AccumPercent AccumDensity AccumDensityOnHash AccumPercent AccumPercentOnHash ISNumInRange FreqInSpecificSize);
use Sort::PureNum;
use BamRelated::BaseAndQual;

# 从单变异对应的插入片段统计表中获得对应的插入片段信息;
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
# 统计（Bed区间、所有区域）内的插入片段分布信息（不分染色体，默认为0~300范围内）；
# >>>>> 这里只会统计数值为正的长度 <<<<<<
sub SizeInfoGetOnBed
{
	my ($Bam,$Samtools,$Bed,$MinSize,$MaxSize) = @_;
	my @InsertInfo = ();
	
	$MinSize = 1 unless($MinSize);
	$MaxSize = 300 unless($MaxSize);
	my $Num = $MaxSize - $MinSize;
	die "[ Error ] MaxSize < MinSize.\n" unless($Num > 0);
	die "[ Error ] Bam or Samtools not exist. ($Bam,$Samtools)\n" unless(-s $Bam && -s $Samtools);
	if($Bed)
	{
		die "[ Error ] Bed not exist. ($Bed)\n" unless(-s $Bed);
	}
	for my $i (0 .. $Num)
	{
		$InsertInfo[$i] = 0;
	}
	
	open(BH,"$Samtools view -F 0xD04 $Bam | cut -f 9 |") or die $! if($Bed);
	open(BH,"$Samtools view -F 0xD04 -L $Bed $Bam | cut -f 9 |") or die $! unless($Bed);
	while(my $Size = <BH>)
	{
		chomp $Size;
		next unless($Size >= $MinSize && $Size <= $MaxSize);
		
		my $Id = $Size - $MinSize;
		$InsertInfo[$Id] ++;
	}
	close BH;
	
	return \@InsertInfo;
}
# 单个区域内的片段分布信息
# >>>>> 需要判断是否会重复统计 <<<<<<<
sub SizeInfoGetOnArea
{
	my ($Bam,$Samtools,$Chr,$From,$To,$MinSize,$MaxSize) = @_;
	my @InsertInfo = ();
	
	$MinSize = 1 unless($MinSize);
	$MaxSize = 300 unless($MaxSize);
	my $Num = $MaxSize - $MinSize;
	die "[ Error ] MaxSize < MinSize.\n" unless($Num > 0);
	die "[ Error ] Bam or Samtools not exist. ($Bam,$Samtools)\n" unless(-s $Bam && -s $Samtools);
	for my $i (0 .. $Num)
	{
		$InsertInfo[$i] = 0;
	}
	
	my $Area = "";
	$Area = "$Chr" . ":" . $From . "-" . $To if($Chr && $From && $To);
	open(BH,"$Samtools view -F 0xD04 $Bam $Area |") or die $! if($Area);
	open(BH,"$Samtools view -F 0xD04 $Bam |") or die $! unless($Area);
	while(my $Line = <BH>)
	{
		chomp $Line;
		my @Cols = split /\t/, $Line;
		my $Size = $Cols[8];
		# 在Size为负时需要判断配对Reads是否会被覆盖到；
		if($Size < 0)
		{
			my $ReadLen = length($Cols[9]);
			# 假如正向会被统计那么负向的忽略就好；
			next if($Cols[7] + $ReadLen > $From);
			$Size = abs($Size);
		}
		next unless($Size >= $MinSize && $Size <= $MaxSize);
		
		my $Id = $Size - $MinSize;
		$InsertInfo[$Id] ++;
	}
	close BH;
	
	for my $i (0 .. $Num)
	{
		$InsertInfo[$i] = 0 unless($InsertInfo[$i]);
	}
	
	return \@InsertInfo;
}
# 统计不同染色体的插入片段分布信息（可以指定bed，默认为0~300范围内）；
# >>>>> 这里只会统计数值为正的长度 <<<<<<
sub SizeInfoGetOnChr
{
	my ($Bam,$Samtools,$Bed,$MinSize,$MaxSize) = @_;
	my %InsertInfo = ();
	
	$MinSize = 1 unless($MinSize);
	$MaxSize = 300 unless($MaxSize);
	my $Num = $MaxSize - $MinSize;
	die "[ Error ] MaxSize < MinSize.\n" unless($Num > 0);
	die "[ Error ] Bam or Samtools not exist. ($Bam,$Samtools)\n" unless(-s $Bam && -s $Samtools);
	if($Bed)
	{
		die "[ Error ] Bed not exist. ($Bed)\n" unless(-s $Bed);
	}
	
	open(BH,"$Samtools view -F 0xD04 -L $Bed $Bam | cut -f 3,9 |") or die $! if($Bed);
	open(BH,"$Samtools view -F 0xD04 $Bam | cut -f 3,9 |") or die $! unless($Bed);
	while(my $Line = <BH>)
	{
		chomp $Line;
		my ($Chr,$Size) = split /\t/, $Line;
		next unless($Size >= $MinSize && $Size <= $MaxSize);
		
		my $Id = $Size - $MinSize;
		$InsertInfo{$Chr}[$Id] = 0 unless($InsertInfo{$Chr}[$Id]);
		$InsertInfo{$Chr}[$Id] ++;
	}
	close BH;
	
	foreach my $Chr (keys %InsertInfo)
	{
		for my $i (0 .. $Num)
		{
			$InsertInfo{$Chr}[$i] = 0 unless($InsertInfo{$Chr}[$i]);
		}
	}
	
	return \%InsertInfo;
}
# 统计给定Bed区间内所有位点的片段分布信息；
# >>>>> 需要判断是否会重复统计 <<<<<<<
sub SizeInfoOfSingleOnBed
{
	my ($Bam,$Samtools,$Bed,$MinSize,$MaxSize) = @_;
	my %SizeInfo = ();
	
	# 在统计单点的条件下，点与点之间是相互独立的，因此不需要考虑重叠的问题;
	$MinSize = 1 unless($MinSize);
	$MaxSize = 300 unless($MaxSize);
	my $Num = $MaxSize - $MinSize;
	die "[ Error ] MaxSize < MinSize.\n" unless($Num > 0);
	die "[ Error ] Bam samtools or bed not specified.\n" unless($Bam && $Samtools && $Bed);
	open(BAM,"$Samtools view -F 0xD04 -L $Bed $Bam |") or die $!;
	while(my $Line = <BAM>)
	{
		chomp $Line;
		my @Cols = split /\t/, $Line;
		next unless(abs($Cols[8]) >= $MinSize && abs($Cols[8]) <= $MaxSize);
		
		# 计算归属时不考虑SoftClip；
		my $SclipFlag = 0;
		my ($Chr,$From,$To) = ReadCoverRange($Line,$SclipFlag);
		# 需要过滤掉正向重叠的区域；
		if($Cols[8] < 0)
		{
			my $ReadLen = length($Cols[9]);
			# 忽略正向会被统计的部分；
			$From = $Cols[7] + $ReadLen if($Cols[7] + $ReadLen > $From);
			$Cols[8] = abs($Cols[8]);
		}
		my $Id = $Cols[8] - $MinSize;
		for my $i ($From .. $To)
		{
			my $Key = join("\t",$Chr,$i);
			$SizeInfo{$Key}[$Id] = 0 unless($SizeInfo{$Key}[$Id]);
			$SizeInfo{$Key}[$Id] ++;
		}
	}
	close BAM;
	
	foreach my $Key (keys %SizeInfo)
	{
		for my $i (0 .. $Num)
		{
			$SizeInfo{$Key}[$i] = 0 unless($SizeInfo{$Key}[$i]);
		}
	}
	
	return \%SizeInfo;
}
sub SizeInfoOfMultiArea
{
	my ($Bam,$Samtools,$SortBed,$MergeBed,$MinSize,$MaxSize) = @_;
	my @SizeInfo = ();
	
	$MinSize = 1 unless($MinSize);
	$MaxSize = 300 unless($MaxSize);
	my $Num = $MaxSize - $MinSize;
	die "[ Error ] MaxSize < MinSize.\n" unless($Num > 0);
	die "[ Error ] Bam samtools or bed not specified.\n" unless($Bam && $Samtools && $SortBed && $MergeBed);
	
	# 获取所有片段信息；
	my (@AChr,@AFrom,@ATo) = ();
	open(BED,"cat $SortBed |") or die $! unless($Bed =~ /\.gz$/);
	open(BED,"zcat $SortBed |") or die $! if($Bed =~ /\.gz$/);
	while(my $Line = <BED>)
	{
		chomp $Line;
		my @Cols = split /\t/, $Line;
		push @AChr, $Cols[0];
		push @AFrom, $Cols[1];
		push @ATo, $Cols[2];
	}
	close BED;
	
	# 从编号0开始顺序搜寻，没有就下一位；
	my $SearchId = 0;
	my $MaxLoopNum = $#AChr;
	open(BAM,"$Samtools view -F 0xD04 -L $MergeBed $Bam |") or die $!;
	while(my $Line = <BAM>)
	{
		chomp $Line;
		my @Cols = split /\t/, $Line;
		next unless(abs($Cols[8]) >= $MinSize && abs($Cols[8]) <= $MaxSize);
		
		# 计算归属时不考虑SoftClip；
		my $SclipFlag = 0;
		my ($Chr,$From,$To) = ReadCoverRange($Line,$SclipFlag);
		# 需要过滤掉正向重叠的区域；
		if($Cols[8] < 0)
		{
			my $ReadLen = length($Cols[9]);
			# 忽略正向会被统计的部分；
			$From = $Cols[7] + $ReadLen if($Cols[7] + $ReadLen > $From);
			next if($From > $To);
			$Cols[8] = abs($Cols[8]);
		}
		my $SizeId = $Cols[8] - $MinSize;
		
		# 定位对应的有重叠的Bed编号；
		my $LoopNum = 0;
		while($AChr[$SearchId] ne $Chr || $AFrom[$SearchId] > $To || $ATo[$SearchId] < $From)
		{
			$SearchId ++;
			$SearchId = 0 if($SearchId > $MaxLoopNum);
			$LoopNum ++;
			last if($LoopNum > $MaxLoopNum);
		}
		
		for my $i ($SearchId .. $#AChr)
		{
			# 有重叠
			if($AChr[$i] eq $Chr && $AFrom[$i] <= $To && $ATo[$i] >= $From)
			{
				$SizeInfo[$i][$SizeId] = 0 unless($SizeInfo[$i][$SizeId]);
				$SizeInfo[$i][$SizeId] ++;
			}
			else
			{
				last;
			}
		}
	}
	close BAM;
	
	for my $i (0 .. $#AChr)
	{
		for my $j (0 .. $Num)
		{
			$SizeInfo[$i][$j] = 0 unless($SizeInfo[$i][$j]);
		}
	}
	
	return \@SizeInfo;
}

# 用于合并固定行号范围内的插入片段数据（文件必须要有注释行且为非压缩文件）；
# 其中Prefix指的就是第一列的信息；
sub SizeCountMerge
{
	my ($File,$Prefix,$IdFrom,$IdTo,$tRef) = @_;
	my @SizeInfo = ();
	@SizeInfo = @{$tRef} if($tRef);
	
	my $InitialId = 0;
	my $Header = "";
	$Header = `cat $File | grep ^# | tail -n1` unless($File =~ /\.gz$/);
	$Header = `zcat $File | grep ^# | tail -n1` if($File =~ /\.gz$/);
	chomp $Header;
	my @Items = split /\t/, $Header;
	$Items[0] =~ s/^#//;
	for my $i (0 .. $#Items)
	{
		unless($Items[$i] =~ /\D/)
		{
			$InitialId = $i;
			last;
		}
	}
	
	$IdFrom = 1 unless($IdFrom);
	unless($IdTo)
	{
		if($File =~ /\.gz$/)
		{
			$IdTo = `zcat $File | grep -v ^# | wc -l` unless($Prefix);
			$IdTo = `zcat $File | grep -v ^# | awk '{if(\$1 == "$Prefix"){print \$0}}' | wc -l` if($Prefix);
		}
		else
		{
			$IdTo = `cat $File | grep -v ^# | wc -l` unless($Prefix);
			$IdTo = `cat $File | grep -v ^# | awk '{if(\$1 == "$Prefix"){print \$0}}' | wc -l` if($Prefix);
		}
		chomp $IdTo;
	}
	die "[ Error ] IdTo < IdFrom.\n" unless($IdTo >= $IdFrom);
	
	die "[ Error ] File not exist ($File).\n" unless(-s $File);
	if($File =~ /\.gz$/)
	{
		open(FI,"zcat $File | grep -v ^# | awk '{if(NR >= $IdFrom && NR <= $IdTo){print \$0}}' |") or die $! unless($Prefix);
		open(FI,"zcat $File | grep -v ^# | awk '{if(\$1 == \"$Prefix\"){print \$0}}' | awk '{if(NR >= $IdFrom && NR <= $IdTo){print \$0}}' |") or die $! if($Prefix);
	}
	else
	{
		open(FI,"cat $File | grep -v ^# | awk '{if(NR >= $IdFrom && NR <= $IdTo){print \$0}}' |") or die $! unless($Prefix);
		open(FI,"cat $File | grep -v ^# | awk '{if(\$1 == \"$Prefix\"){print \$0}}' | awk '{if(NR >= $IdFrom && NR <= $IdTo){print \$0}}' |") or die $! if($Prefix);
	}
	while(my $Line = <FI>)
	{
		chomp $Line;
		my @Size = split /\t/, $Line;
		for my $i ($InitialId .. $#Size)
		{
			$SizeInfo[$i - $InitialId] = 0 unless($SizeInfo[$i - $InitialId]);
			$SizeInfo[$i - $InitialId] += $Size[$i];
		}
	}
	close FI;
	
	return \@SizeInfo;
}
# 用于根据Lo et al., 2014 PNAS中的方式计算FF；
sub SizeCountBasedOnLo2014
{
	my $File = $_[0];
	my @SizeInfo = @{$_[1]};
	
	# 确定数字编码起始；
	my $InitialId = -1;
	my $Header = "";
	$Header = `cat $File | grep ^# | tail -n1` unless($File =~ /\.gz$/);
	$Header = `zcat $File | grep ^# | tail -n1` if($File =~ /\.gz$/);
	chomp $Header;
	my @HeaderItems = split /\t/, $Header;
	$HeaderItems[0] =~ s/^#//;
	for my $i (0 .. $#HeaderItems)
	{
		unless($HeaderItems[$i] =~ /\D/)
		{
			$InitialId = $i;
			last;
		}
	}
	
	# 统计100到150；
	my $NumOfShort = 0;
	for my $i (100 .. 150)
	{
		my $ColId4Num = 0;
		for my $j ($InitialId .. $#HeaderItems)
		{
			if($HeaderItems[$j] == $i)
			{
				$ColId4Num = $j;
				last;
			}
		}
		die "[ Error ] Cannot locate $j in the header of $File\.\n" if($ColId4Num < $InitialId);
		
		$NumOfShort += $SizeInfo[$ColId4Num - $InitialId];
	}
	
	# 统计163到169；
	my $NumOfLong = 0;
	for my $i (163 .. 169)
	{
		my $ColId4Num = 0;
		for my $j ($InitialId .. $#HeaderItems)
		{
			if($HeaderItems[$j] == $i)
			{
				$ColId4Num = $j;
				last;
			}
		}
		die "[ Error ] Cannot locate $j in the header of $File\.\n" if($ColId4Num < $InitialId);
		
		$NumOfLong += $SizeInfo[$ColId4Num - $InitialId];
	}
	
	die "[ Error ] Num from 163 to 169 was zero (",join("\t",@SizeInfo),")\n" unless($NumOfLong > 0);
	my $FF = $NumOfShort / $NumOfLong;
	
	return ($FF,$NumOfShort,$NumOfLong);
}
# 用于两个分布相减并给出最大值；
sub SubtractOfDistrMax
{
	my @AltDistr = @{$_[0]};
	my @RefDistr = @{$_[1]};
	
	for my $i (0 .. $#AltDistr)
	{
		$AltDistr[$i] = $AltDistr[$i] - $RefDistr[$i];
	}
	
	my $Max = 0;
	for my $i (0 .. $#AltDistr)
	{
		$Max = $AltDistr[$i] if(abs($AltDistr[$i]) > abs($Max));
	}
	
	return $Max;
}
# 用于两个分布相减并给出平均值；
sub SubtractOfDistrMean
{
	my @AltDistr = @{$_[0]};
	my @RefDistr = @{$_[1]};
	
	for my $i (0 .. $#AltDistr)
	{
		$AltDistr[$i] = $AltDistr[$i] - $RefDistr[$i];
	}
	
	my $Sum = 0;
	my $Num = 0;
	for my $i (0 .. $#AltDistr)
	{
		$Sum += $AltDistr[$i];
		$Num ++;
	}
	my $Mean = $Sum / $Num;
	
	return $Mean;
}

# 将数量信息转换成比例信息；
sub Size2Percent
{
	my @Size = @{$_[0]};
	my @Percent = ();
	
	my $Total = 0;
	for my $i (0 .. $#Size)
	{
		$Total += $Size[$i];
	}
	for my $i (0 .. $#Size)
	{
		$Percent[$i] = $Size[$i] / $Total;
	}
	
	return \@Percent;
}
sub Size2AccumPercent
{
	my @Size = @{$_[0]};
	my @Percent = ();
	
	my $Total = 0;
	for my $i (0 .. $#Size)
	{
		$Total += $Size[$i];
	}
	my $Num = 0;
	for my $i (0 .. $#Size)
	{
		$Num += $Size[$i];
		$Percent[$i] = $Num / $Total;
	}
	
	return \@Percent;
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
