#!/usr/bin/perl
use Getopt::Long;
&GetOptions('h=s' =>\$helpFlag,
'fai=s' =>\$inputFile1,
'file=s' =>\$inputFile2,
'type=s' =>\$type);
$usage="Correct-cufflinks-predition-for-chr-length-V1.pl 
-fai fai file
-file gtf file name (bed/gtf/gff)
-type bed,gtf,gff
nimarafati\@gmail.com 20150323\n";
if($helpFlag || $inputFile1 eq ""|| $inputFile2 eq "")
{
	print $usage;
	exit;
}
open(inF1,$inputFile1);
while(<inF1>)
{
#	print $_;<STDIN>;
	chomp($_);
	@lineArr=split("\t",$_);
	$chrHash{$lineArr[0]}=$lineArr[1];
#	print "here:",$lineArr[$col_key-1],"\t$chrHash{$lineArr[$col_key-1]}";<STDIN>;
}
close inF1;
open(inF2,$inputFile2);
while(<inF2>)
{
	chomp($_);
	@lineArr1=split("\t",$_);
	if($type eq "gtf" || $type eq "gff")
	{
		$end=$lineArr1[4];
	}
	else
	{
		$end=$lineArr1[2];
	}
	if(exists $chrHash{$lineArr1[0]})
	{
		for(my $i=0;$i<scalar(@lineArr1);$i++)
		{
			if($lineArr1[$i]>$chrHash{$lineArr1[0]})
			{
				print $chrHash{$lineArr1[0]},"\t";
			}
			else
			{
				print $lineArr1[$i],"\t";
			}
		}
		print "\n";
	}
}
close inF2;
