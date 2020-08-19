#!/usr/bin/perl
#
use Getopt::Long;

my $contact="Nima Rafati nimarafati\@gmail.com
2016	";
my $usage = "./*.pl -i file\n$contact\n";

&GetOptions('h=s' =>\$helpFlag,
'i=s' =>\$inputFile);
if($helpFlag)
{
	print $usage;
	exit;
}

my @lineArr=();
my $cntr=0;
open(inF1,$inputFile) || die print $usage;
while(<inF1>)
{
	$cntr++;
	$cntr_1=0;
	chomp($_);
	@lineArr=split("\t",$_);
	for($i=3;$i<scalar(@lineArr);$i++)
	{
		if($lineArr[$i]==1)
		{
			$cntr_1++;
		}
		if($cntr_1==3)
		{
			print $_,"\n";
		}
	}
}
close inF1;

