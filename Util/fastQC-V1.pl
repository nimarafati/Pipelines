#!/usr/bin/perl
#This script gets paths and generate fastQC commands

use Getopt::Long;
use File::Basename;
my $contact="Nima Rafati nimarafati\@gmail.com
20141104	";
my $usage="\nfastQC-V1.pl -path-files file.txt -t 2
-path-files a file including path to directory where there are fastq files
-t number of threads max = 6
$contact\n";
$threads=2;
&GetOptions('path-files=s' =>\$listFile, 'h=s' =>\$helpFlag , 't=i' =>\$threads);

if($listFile eq "" || $helpFlag)
{
	print $usage;
	exit;
}
if($threads >6)
{
	$threads=6;
}
open(inF1,$listFile);
while(<inF1>)
{
	chomp($_);
	$line=$_;
	my($file, $dir, $ext) = fileparse($line);
	opendir($dh,$line);
	my @list=readdir($dh);
#	my @list= grep { /^\./ && -f "$line" } readdir($dh);
	print "mkdir $file\nfastqc --contaminants ~/glob/Contamination/adapters-tab-txt.list -o $file --noextract -t $threads ";
	foreach my $i (@list)
	{
		if($i=~ m/^\.\./ )#|| $i ne "..")
		{
#			print "1";
		}
		elsif($i=~ m/^\./)
		{
#			print "2";
		}
		else
		{
			print " $dir$file/$i";
		}
	}
	print "\n";
	closedir($dh);
}
close inF1;
