#!/usr/bin/perl
open(inF1,$ARGV[0]);
while(<inF1>)
{
	chomp($_);
	@lineArr=split("\t",$_);
	$chrHash{$lineArr[0]}=$lineArr[0];
}
close inF1;
open(inF2,$ARGV[1]);
while(<inF2>)
{
	chomp($_);
	@lineArr1=split("\t",$_);
	if(exists $chrHash{$lineArr1[0]})
	{
		for ($i=0;$i<scalar(@lineArr1);$i++)
		{
			print "$lineArr1[$i]\t";
		}
		print "\n";
	}
}
close inF2;

