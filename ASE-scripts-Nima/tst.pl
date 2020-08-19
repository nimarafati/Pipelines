#!/usr/bin/perl

#$par1= system "echo test1";
#$par2= system "echo test2";
open(inF1,"samtools view ../Sample_AB-d012-AA/Sample_AB-d012-gsnap-snpT-AA.Reordered.sort.bam |");
open(inF2,"samtools view ../Sample_AB-d012-AA/Sample_AB-d012-gsnap-snpT-AA.Reordered.sort.bam |");
	while($par1 = <inF1>)
	{
		my $par1=<inF1>;
		my $par2=<inF2>;
		print "1:$par1\n2:$par2";<STDIN>;
		@lineArr=split("\t",$par1);
		$t=cigar(@lineArr);
	}
close inF1;
close inF2;
sub cigar {
        print join "\t",@_ if $verbose;
        # This should be improved since the column order may differ. However, the current implementation works with GSNAP.
foreach my $ar (@_){if ($ar=~ m/MD.*/){$MD=$ar;}}
my $cigar=$_[5];my $qual=$_[10];
print "$MD\n$qual\n$cigar";<STDIN>;

}
