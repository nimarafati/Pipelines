#!/usr/bin/perl
use Cwd;
use File::Find;
my $cWD=getcwd;
my $cntr=0;
open(outR,">Samples-QC-summary");
open(inSample,"Samples.list");
while(<inSample>)
{
        chomp($_);
        $line_sample=$_;
        @lineArr=split("\t",$_);
        chdir "$cWD/fastQC/$lineArr[0]";
#	print getcwd;<STDIN>;
	opendir(DIR,"$cWD/fastQC/$lineArr[0]");
	my @files = grep {!/^\./} readdir DIR;
	foreach my $f (@files)
	{
		$f=~ m/(.*)\.zip/;
		$zip_file=$1;
		if($zip_file ne "" && $cntr==0)
		{
			chdir "$cWD/fastQC/$lineArr[0]/$zip_file";	
			system "unzip -fq $f";
	        	open (inFastQC_data,"fastqc_data.txt");
		        while(<inFastQC_data>)
	        	{
	        	        chomp($_);
        		        $line_qc1=$_;
	        	        @lineArr_qc1=split("\t",$_);
        		        foreach my $q (@lineArr_qc1)
	                	{
	                	        if($q eq "fail")
        		                {
        	        	                $lineArr_qc1[0]=~ m/>>(.*)/;
	                        	        $failed1.="$1,";
	                        	}
        	        	}
				if($line_qc1=~  m/Total Sequences.*/){
				$line_qc1=~ m/(Total)\s(Sequences)\s(\d+)/;
		                $read_count1=$3;
			#	print "RC:$read_count";<STDIN>;
				}
	        	}
			close inFastQC_data;
			$cntr=1;
		}
	}
	$cntr=0;
        $report.="$line_sample\t$read_count1\t$failed1\t";
        $failed1="";
        $read_count1="";
	@files=();
	chdir "$cWD/fastQC-after-trimming/$lineArr[0]";
	opendir(DIR,"$cWD/fastQC-after-trimming/$lineArr[0]");
	my @files = grep {!/^\./} readdir DIR;
#	print @files;<STDIN>;
	foreach my $f (@files)
	{
		$f=~ m/(.*)\.zip/;
		$zip_file=$1;
		if($zip_file ne "" && $cntr==0)
		{
			chdir "$cWD/fastQC-after-trimming/$lineArr[0]/$zip_file";
#			print getcwd",$f";<STDIN>;
#			print "UNZIP this one:\n$f in ",getcwd ;<STDIN>;
			system "unzip -q $f";
	        	open (inFastQC_data,"fastqc_data.txt");
		        while(<inFastQC_data>)
	        	{
	        	        chomp($_);
        		        $line_qc2=$_;
	        	        @lineArr_qc2=split("\t",$_);
         		        foreach my $q (@lineArr_qc2)
	                	{
	                	        if($q eq "fail")
        		                {
        	        	                $lineArr_qc2[0]=~ m/>>(.*)/;
	                        	        $failed2.="$1,";
	                        	}
        	        	}
				if($line_qc2=~  m/Total Sequences.*/){
				$line_qc2=~ m/(Total)\s(Sequences)\s(\d+)/;
		                $read_count2=$3;
			#	print "RC:$read_count";<STDIN>;
				}
	        	}
			close inFastQC_data;
			$cntr=1;
		}
	}
	$cntr=0;
	$report.="$read_count2\t$failed2\t";
        $failed2="";
        $read_count2="";
	print outR "$report\n";
	$report="";
	$zip_file="";
}
close inSample;
close outR;
