#!/usr/bin/perl
#This script generates commands for trimmomatic

use Getopt::Long;
use File::Basename;
use Cwd;
my $trimPath="~/glob/Software/Trimmomatic-0.30/trimmomatic-0.30.jar";
my $adapterFile="/home/nimar/glob/Contamination/adapters.fa";
my $qualCut=15;
my $minLength=36;
my $phred=33;
my $cwd=getcwd;
my $contact="Nima Rafati nimarafati\@gmail.com
20141104	";
my $usage="\nfastQC-V1.pl -path-files file.txt -t 2 -trimmomatic-path -adapter-file -q -MINLEN -phred
-path-files a file including path to directory where there are fastq files
-t number of threads max = 6
-trimmomatic-path path to trimmomatic (default=~/glob/Software/Trimmomatic-0.30/trimmomatic-0.30.jar)
-adapter-file adapter.fa file (default=/home/nimar/glob/Contamination/adapters.fa)
-q quality threshols (default=15)
-MINLEN minimum length (default=36)
-phred phred scale (default=33)
$contact\n";
$threads=2;
&GetOptions('path-files=s' =>\$listFile, 'h=s' =>\$helpFlag , 't=i' =>\$threads, 'trimmomatic-path=s' =>\$trimPath, 'adapter-file=s' =>\$adapterFile, 'q=i' =>\$qualityCut, 'MINLEN=i' =>\$minLength, 'phred=i' =>\$phred);

if($listFile eq "" || $helpFlag || $adapterFile eq "")
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
	print "mkdir $cwd/$file
cd $cwd/$file
cp $adapterFile adapters.fa
java -jar $trimPath PE -threads $threads -phred$phred -trimlog $file.log ";
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
##java -jar ~/glob/Software/Trimmomatic-0.30/trimmomatic-0.30.jar PE -threads 4 -phred33 -trimlog Sample_48B-trimm.log    /home/nimar/b2013097/private/UCDxJFtranscriptome/130822_D00118_0112_BD2D4AACXX/Sample_48B/48B_ATCACG_L004_R1_001.fastq.gz /home/nimar/b2013097/private/UCDxJFtranscriptome/130822_D00118_0112_BD2D4AACXX/Sample_48B/48B_ATCACG_L004_R2_001.fastq.gz 48B_ATCACG_L004_R1_001-paired.fastq 48B_ATCACG_L004_R1_001-unpaire.fastq 48B_ATCACG_L004_R2_001-paired.fastq 48B_ATCACG_L004_R2_001-unpaired.fastq ILLUMINACLIP:adapters.fa:2:40:15:8:true SLIDINGWINDOW:4:15 LEADING:15 TRAILING:15 MINLEN:36
			print " $dir$file/$i ";
		}
	}
	foreach my $i (@list)
	{
		if($i=~ m/^\.\./ )#|| $i ne "..")
		{
		}
		elsif($i=~ m/^\./)
		{
		}
		else
		{
			$i=~ s/\.fastq.gz//;
			print " $i-paired.fq $i-unpaired.fq ";
		}
	}
	print "ILLUMINACLIP:adapters.fa:2:40:15:8:true SLIDINGWINDOW:4:15 LEADING:15 TRAILING:15 MINLEN:36 &\n";
	closedir($dh);
}
close inF1;
