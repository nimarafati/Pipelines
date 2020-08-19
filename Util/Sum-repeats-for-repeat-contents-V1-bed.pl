#!/usr/bin/perl
#This script sums the repeats' proportion for each annotation and reports that in last column of bed file.
$prv_gene_id="";
open(inF1,$ARGV[0]) || die print $usage;
while(<inF1>)
{
	$cntr++;
	chomp($_);
	@lineArr=split("\t",$_);
	if ($cntr>=1 && $lineArr[0] ne "CHROM")
	{
		chomp($_);
		$gene_id=$lineArr[3];
		$size=$lineArr[2]-$lineArr[1];
		$repeatSize=$lineArr[scalar(@lineArr)-1];
		if($prv_gene_id eq $gene_id || $prv_gene_id eq "")# || $cntr==1)
		{	
			$sumRepeatRatio+=sprintf("%.3f",$repeatSize/$size);
			print "$_,\n$gene_id\t$prv_gene_id\t$repeatSize/$size\t$sumRepeatRatio";<STDIN>;
		}
		else
		{
			@prvLineArr=split("\t",$prv_line);
			for(my $i=0;$i<12;$i++)
			{
				print "$prvLineArr[$i]\t";
			}
			print $sumRepeatRatio,"\n";
			$sumRepeatRatio=0;
			$sumRepeatRatio+=sprintf("%.3f",$repeatSize/$size);
		}
	}
	$prv_line=$_;
	$prv_gene_id=$gene_id;
}
close inF1;
                        @prvLineArr=split("\t",$prv_line);
                        for(my $i=0;$i<12;$i++)
                        {
                                print "$prvLineArr[$i]\t";
                        }
                        print $sumRepeatRatio,"\n";
                        $sumRepeatRatio=0;
                        $sumRepeatRatio+=sprintf("%.3f",$repeatSize/$size);
               

