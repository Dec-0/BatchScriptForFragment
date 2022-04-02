# Package Name
package BedRelated::Bed;

# Exported name
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(BedNarrowDown BedSize SNPListOnBed);

use FindBin qw($Bin);
use File::Basename;

# 用于随机抽取一部分点，生成新的bed文件;
sub BedNarrowDown
{
	my ($BedGen,$OriBed,$RandBed,$RandNum) = @_;
	
	my $Total = 0;
	$Total = `cat $OriBed | awk 'BEGIN{SUM = 0}{SUM += \$3 - \$2}END{print SUM}'` unless($OriBed =~ /\.gz$/);
	$Total = `zcat $OriBed | awk 'BEGIN{SUM = 0}{SUM += \$3 - \$2}END{print SUM}'` if($OriBed =~ /\.gz$/);
	chomp $Total;
	die "[ Error ] Not enough points in $OriBed\n" if($Total < $RandNum);
	
	# 随机抽取;
	my %Points = ();
	if(1)
	{
		# 获得基本信息;
		my @BedItem = ();
		my $LineId = 0;
		open(BED,"cat $OriBed |") or die $! unless($OriBed =~ /\.gz$/);
		open(BED,"zcat $OriBed |") or die $! if($OriBed =~ /\.gz$/);
		while(my $Line = <BED>)
		{
			chomp $Line;
			my @Cols = split /\t/, $Line;
			
			$BedItem[$LineId][0] = $Cols[0];
			$BedItem[$LineId][1] = $Cols[1];
			$BedItem[$LineId][2] = $Cols[2] - $Cols[1];
			$LineId ++;
		}
		close BED;
		
		# 挑选位点;
		for my $i (1 .. $RandNum)
		{
			my $Key = "";
			do
			{
				my $LId = int(rand($LineId));
				my $Pos = $BedItem[$LId][1] + int(rand($BedItem[$LId][2]));
				$Key = join("\t",$BedItem[$LId][0],$Pos);
			}
			while($Points{$Key});
			
			$Points{$Key} = 1;
		}
	}
	
	my $RandFile = $RandBed;
	$RandFile =~ s/bed$/xls/;
	open(TMP,"> $RandFile") or die $!;
	foreach my $Key (keys %Points)
	{
		print TMP $Key,"\n";
	}
	close TMP;
	
	`$BedGen -i $RandFile -o $RandBed`;
	`rm $RandFile`;
	
	return 1;
}

# 计算Bed文件大小;
sub BedSize
{
	my $Bed = $_[0];
	
	die "[ Error ] Bed ($Bed) not exist or empty.\n" unless(-s $Bed);
	my $BedSize = 0;
	$BedSize = `cat $Bed | awk 'BEGIN{Sum = 0}{Sum += \$3 - \$2}END{print Sum}'` unless($Bed =~ /\.gz$/);
	$BedSize = `zcat $Bed | awk 'BEGIN{Sum = 0}{Sum += \$3 - \$2}END{print Sum}'` if($Bed =~ /\.gz$/);
	chomp $BedSize;
	
	return $BedSize;
}

# 返回bed文件对应的snp位点列表文件；
sub SNPListOnBed
{
	my ($Bed) = @_;
	
	my $SNPList = "";
	# 假如该Bed文件就是Bed列表，就不需要继续搜寻；
	$SNPList = $Bed if($Bed =~ /\.MafAll\.bed$/);
	
	# 检测同一目录下是否存在对应的snp列表文件；
	unless($SNPList)
	{
		my $tBed = $Bed;
		$tBed =~ s/bed$/MafAll\.bed/;
		$SNPList = $tBed if(-s $tBed);
	}
	
	return $SNPList;
}

1;