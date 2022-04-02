# Package Name
package MGDRelated::InsertSize2;

# Exported name
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(SizeInfoOfSingleOnBed2 SizeInfoOfMultiArea2);
use BamRelated::BaseAndQual;

# 统计给定Bed区间内所有单个位点对应的片段信息，由于有了SizeInfoOfSingleOnBed因此名称为SizeInfoOfSingleOnBed2。
sub SizeInfoOfSingleOnBed2
{
	my ($Bam,$Samtools,$Bed,$MinSize,$MaxSize) = @_;
	
	# bam、samtools和bed文件需要指定；
	die "[ Error ] Bam samtools or bed not specified.\n" unless($Bam && $Samtools && $Bed);
	
	# 需要纳入统计的片段长度范围；
	$MinSize = 1 unless($MinSize);
	$MaxSize = 300 unless($MaxSize);
	my $SizeRange = $MaxSize - $MinSize;
	die "[ Error ] MinSize ($MinSize) <= 0.\n" unless($MinSize > 0);
	die "[ Error ] MaxSize < MinSize.\n" unless($SizeRange > 0);
	
	# bed文件中同一条染色体对应的坐标需要经过由小到大排序，而染色体号对应的顺序可以不排序。
	
	# 记录Bed文件相关位点，Flag4PosPos中前一个键的值记录后一个的键，Flag4PrePos中后一个键的值记录前一个的键；
	my %Flag4Pos = ();
	my $PreKey = "";
	open(BED,"cat $Bed | cut -f 1-3 |") or die $! unless($Bed =~ /\.gz$/);
	open(BED,"zcat $Bed | cut -f 1-3 |") or die $! if($Bed =~ /\.gz$/);
	while(my $Line = <BED>)
	{
		chomp $Line;
		my @Cols = split /\t/, $Line;
		for my $i ($Cols[1] + 1 .. $Cols[2])
		{
			my $Key = join("\t",$Cols[0],$i);
			$Flag4Pos{$Key} = "None";
			$Flag4Pos{$PreKey} = $Key if($PreKey);
			$PreKey = $Key;
		}
	}
	close BED;
	
	# 逐行处理bed区域内相关reads；
	my %SizeInfo = ();
	open(BAM,"$Samtools view -F 0xD04 -L $Bed $Bam | cut -f 1-10 |") or die $!;
	while(my $Line = <BAM>)
	{
		chomp $Line;
		my @Cols = split /\t/, $Line;
		next unless(abs($Cols[8]) >= $MinSize && abs($Cols[8]) <= $MaxSize);
		
		# 检查该Read对应的覆盖范围（SoftClip不纳入覆盖范围）；
		my $SclipFlag = 0;
		my ($Chr,$From,$To) = ReadCoverRange($Line,$SclipFlag);
		# 片段长度为正时可以直接纳入统计范围。
		# 片段长度为负时需要考虑正向Reads是否已经或即将被统计，避免重复统计。
		# 当片段长度短于Reads读长时，由于SoftClip的作用（匹配错位），逆序会在顺序之前出现。
		my ($SegFrom,$SegTo) = ($From,$To);
		if($Cols[8] > 0)
		{
			$SegTo = $SegFrom + $Cols[8] - 1;
		}
		else
		{
			# 判断正向序列是否已经被统计过或者即将被统计，也就是是否在相关位点列表内；
			my $PairFrom = $Cols[7];
			my $PairTo = $PairFrom - $Cols[8] - 1;
			my $MatchFlag = 0;
			for my $i ($PairFrom .. $PairTo)
			{
				my $Key = join("\t",$Chr,$i);
				if($Flag4Pos{$Key})
				{
					$MatchFlag = 1;
					last;
				}
			}
			next if($MatchFlag);
			
			$SegFrom = $SegTo + $Cols[8] + 1;
		}
		
		# 需要记录的下标；
		my $Id = $SegTo - $SegFrom + 1 - $MinSize;
		
		# 找到对应的第一个点
		my $CurrentKey = "";
		for my $i ($SegFrom .. $SegTo)
		{
			my $Key = join("\t",$Chr,$i);
			if($Flag4Pos{$Key})
			{
				$CurrentKey = $Key;
				last;
			}
		}
		# 跳跃定位；
		while($CurrentKey)
		{
			$SizeInfo{$CurrentKey}[$Id] = 0 unless($SizeInfo{$CurrentKey}[$Id]);
			$SizeInfo{$CurrentKey}[$Id] ++;
			
			$CurrentKey = $Flag4Pos{$CurrentKey};
			last if($CurrentKey eq "None");
			my ($tChr,$tTo) = split /\t/, $CurrentKey;
			last if($tChr ne $Chr || $tTo > $SegTo);
		}
	}
	close BAM;
	
	# 没有被统计到的点也需要记录为0;
	foreach my $Key (keys %Flag4Pos)
	{
		for my $i (0 .. $SizeRange)
		{
			$SizeInfo{$Key}[$i] = 0 unless($SizeInfo{$Key}[$i]);
		}
	}
	
	return \%SizeInfo;
}

# 统计给定Bed区间内指定区域对应的片段信息，由于有了SizeInfoOfMultiArea因此名称为SizeInfoOfMultiArea2。
sub SizeInfoOfMultiArea2
{
	# File4Mark记录Bed文件中每行bed对应的标签，假如不指定该文件，那么以染色体号为准；
	my ($Bam,$Samtools,$Bed,$File4Mark,$MinSize,$MaxSize) = @_;
	
	# bam、samtools和bed文件需要指定；
	die "[ Error ] Bam samtools or bed not specified.\n" unless($Bam && $Samtools && $Bed);
	# 用于记录Bed行对应对应统计单位的文件；
	if($File4Mark)
	{
		die "[ Error ] File not exist ($File4Mark).\n" unless(-s $File4Mark);
	}
	
	# 需要纳入统计的片段长度范围；
	$MinSize = 1 unless($MinSize);
	$MaxSize = 300 unless($MaxSize);
	my $SizeRange = $MaxSize - $MinSize;
	die "[ Error ] MinSize ($MinSize) <= 0.\n" unless($MinSize > 0);
	die "[ Error ] MaxSize < MinSize.\n" unless($SizeRange > 0);
	
	# bed文件中同一条染色体对应的坐标需要经过由小到大排序，而染色体号对应的顺序可以不排序。
	
	# 记录Bed文件相关位点，4个数组分别记录chr、from、to和mark，1个哈希记录基因坐标对应的数组下标，1个哈希记录对应的区域标记。
	my (%Flag4Pos,%MarkInfo) = ();
	my (@BChr,@BFrom,@BTo,@BMark) = ();
	open(BED,"cat $Bed | cut -f 1-3 |") or die $! unless($Bed =~ /\.gz$/);
	open(BED,"zcat $Bed | cut -f 1-3 |") or die $! if($Bed =~ /\.gz$/);
	open(MARK,"cat $File4Mark | cut -f 4 |") or die $! if($File4Mark && $File4Mark !~ /\.gz$/);
	open(MARK,"zcat $File4Mark | cut -f 4 |") or die $! if($File4Mark && $File4Mark =~ /\.gz$/);
	while(my $Line = <BED>)
	{
		chomp $Line;
		my @Cols = split /\t/, $Line;
		$Cols[1] ++;
		# 默认情况下标记为染色体号；
		my $Mark = $Cols[0];
		if($File4Mark)
		{
			$Mark = <MARK>;
			chomp $Mark;
		}
		push @BChr, $Cols[0];
		push @BFrom, $Cols[1];
		push @BTo, $Cols[2];
		push @BMark, $Mark;
		$MarkInfo{$Mark} = 1;
		
		# 记录所属数组下标；
		my $tId = @BChr;
		for my $i ($Cols[1] .. $Cols[2])
		{
			my $Key = join("\t",$Cols[0],$i);
			$Flag4Pos{$Key} = $tId;
		}
	}
	close BED;
	close MARK if($File4Mark);
	
	# 逐行处理bed区域内相关reads；
	my %SizeInfo = ();
	open(BAM,"$Samtools view -F 0xD04 -L $Bed $Bam | cut -f 1-10 |") or die $!;
	while(my $Line = <BAM>)
	{
		chomp $Line;
		my @Cols = split /\t/, $Line;
		next unless(abs($Cols[8]) >= $MinSize && abs($Cols[8]) <= $MaxSize);
		
		# 检查该Read对应的覆盖范围（SoftClip不纳入覆盖范围）；
		my $SclipFlag = 0;
		my ($Chr,$From,$To) = ReadCoverRange($Line,$SclipFlag);
		# 片段长度为正时可以直接纳入统计范围。
		# 片段长度为负时需要考虑正向Reads是否已经或即将被统计，避免重复统计。
		# 当片段长度短于Reads读长时，由于SoftClip的作用（匹配错位），逆序会在顺序之前出现。
		my ($SegFrom,$SegTo) = ($From,$To);
		if($Cols[8] > 0)
		{
			$SegTo = $SegFrom + $Cols[8] - 1;
		}
		else
		{
			# 判断正向序列是否已经被统计过或者即将被统计，也就是是否在相关位点列表内；
			my $PairFrom = $Cols[7];
			my $PairTo = $PairFrom - $Cols[8] - 1;
			my $MatchFlag = 0;
			for my $i ($PairFrom .. $PairTo)
			{
				my $Key = join("\t",$Chr,$i);
				if($Flag4Pos{$Key})
				{
					$MatchFlag = 1;
					last;
				}
			}
			next if($MatchFlag);
			
			$SegFrom = $SegTo + $Cols[8] + 1;
		}
		
		# 需要记录的下标；
		my $Id = $SegTo - $SegFrom + 1 - $MinSize;
		
		# 获得第一个数组下标;
		my $BId = 0;
		for my $i ($SegFrom .. $SegTo)
		{
			my $Key = join("\t",$Chr,$i);
			if($Flag4Pos{$Key})
			{
				$BId = $Flag4Pos{$Key};
				last;
			}
		}
		# 顺序检查；
		if($BId > 0)
		{
			# 以防重复标记；
			my %Flag4Mark = ();
			for my $i ($BId - 1 .. $#BChr)
			{
				next if($Flag4Mark{$BMark[$i]});
				last if($BChr[$i] ne $Chr || $SegTo < $BFrom[$i]);
				
				$SizeInfo{$BMark[$i]}[$Id] = 0 unless($SizeInfo{$BMark[$i]}[$Id]);
				$SizeInfo{$BMark[$i]}[$Id] ++;
				$Flag4Mark{$BMark[$i]} = 1;
			}
		}
	}
	close BAM;
	
	# 没有被统计到的点也需要记录为0;
	foreach my $Mark (keys %MarkInfo)
	{
		for my $i (0 .. $SizeRange)
		{
			$SizeInfo{$Mark}[$i] = 0 unless($SizeInfo{$Mark}[$i]);
		}
	}
	
	return \%SizeInfo;
}

1;
