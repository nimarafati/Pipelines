#!/usr/bin/perl
#This script prepares directories to setup RNA-seq pipeline
use Cwd;
use File::Find;
use Getopt::Long;
use POSIX qw(strftime);
$d=system 'date';
$d=strftime "%Y%m%d-%H-%M", localtime();

my $contact="Nima Rafati nimarafati\@gmail.com 
20141012 V1
20150905 V2";

my $usage="Prepare-samples-V3-SE-PE-flag.pl -path 
-path Path to sample files
-REF Genome assembly
-GTF annotation
-VCF VCF file (optional)
-proj project name (Optional)
-readLength length of the reads
-data_type SE/PE single-end paired-end
To update:
1- fastQC report
2- QC after alignment: QoRTs and TIN
$contact";
my $cWD=getcwd;
&GetOptions(
'path=s' =>\$path,
'REF=s' =>\$reference,
'GTF=s' =>\$gtf,
'VCF=s' =>\$vcf,
'proj=s' =>\$proj,
'readLength=s' =>\$readLength,
'data_type=s' => \$data_type);
#'reference_index=s' => \$reference_fai);

if($path eq "" || $reference eq "" || $readLength eq "" || $data_type eq "")
{
	print "Please privde all needed information \n$usage\n";
	exit;
}
system 'mkdir -p fastQC';
system 'mkdir -p Trimmed-reads';
system 'mkdir -p fastQC-after-trimming';
system 'mkdir -p GSNAP';
opendir(DIR, $path);
#print "$path";<STDIN>;
my @files = grep {!/^\./} readdir DIR;
open(outList,">Samples.list");
foreach $dir (sort @files)
{
	
	$dir=~ m/(Sample_.*)/;
	$now=$1;
#	print "DIR=$dir\n$now";<STDIN>;
	if($now ne "" && $now=~ m/Sample_.*/)
	{
		$sampleDir=$path."$now/";
		$fastQC="$cWD/fastQC/$now";
		$Trimmed="$cWD/Trimmed-reads/$now";
		$fastQC_After="$cWD/fastQC-after-trimming/$now";
		$GSNAP="$cWD/GSNAP/$now";
		no_Print_CMD("mkdir -p $GSNAP");
		no_Print_CMD("mkdir -p $fastQC");
		no_Print_CMD("mkdir -p $Trimmed");
		no_Print_CMD("mkdir -p $fastQC_After");
		no_Print_CMD("mkdir -p $GSNAP");
		no_Print_CMD("mkdir -p $cWD/cufflinks");
		opendir(DIR1,$sampleDir);
		my @reads= grep {!/^\./} readdir DIR1;
		foreach $r (sort @reads)
		{
			if(scalar(@reads)>1)
			{
#			print "$sampleDir $r";<STDIN>;
				$readsFile.="\t$r";
			}
			elsif(scalar(@reads==1))
			{
				$readsFile="\t$r\t";
			}
		}
		print outList "$now$readsFile\t$sampleDir\n";
		$readsFile="";
		$sampleDir="";
		$now="";
	}
}
close outList;
##Quality check by FastQC
my $QC1="module load FastQC\n";
$QC1.="cd $cWD/fastQC/\n";
open(inSample,"Samples.list");
while(<inSample>)
{
	$cntr++;
	chomp($_);
	@lineArr=split("\t",$_);
	if($data_type eq "PE")
	{
		$QC1.="fastqc --extract --contaminants /proj/b2012111/private/UserDirectory/nima/glob_old/Contamination/adapters-tab-txt.list -o $lineArr[0]  -t 2 $lineArr[3]/$lineArr[1] $lineArr[3]/$lineArr[2] &\n";
	}
	else
	{
		$QC1.="fastqc --extract --contaminants /proj/b2012111/private/UserDirectory/nima/glob_old/Contamination/adapters-tab-txt.list -o $lineArr[0]  -t 2 $lineArr[3]/$lineArr[1] &\n";
	}
}
Sbatch_run($QC1,"QC-1",$proj,$cntr*40);
$QC1="";
close inSample;
##Trimming by trimmomatic
open(inSample,"Samples.list");
while(<inSample>)
{
	chomp($_);
	@lineArr=split("\t",$_);
	$trim.="cd $cWD/Trimmed-reads/$lineArr[0]\n";
	$trim.="cp /proj/b2012111/private/UserDirectory/nima/glob_old/Contamination/adapters.fa adapters.fa \n";
	$lineArr[1]=~ m/(.*)\.fastq.gz/;
	$for=$1;
	if($data_type eq "PE")
	{
		$lineArr[2]=~ m/(.*)\.fastq.gz/;
		$rev=$1;
		$trim.="java -jar /proj/b2012111/private/UserDirectory/nima/glob_old/Software/Trimmomatic-0.30/trimmomatic-0.30.jar PE -threads 2 -phred33 -trimlog $lineArr[0].log  $lineArr[3]/$lineArr[1] $lineArr[3]/$lineArr[2] $for-paired.fq $for-unpaired.fq $rev-paired.fq $rev-unpaired.fq ILLUMINACLIP:adapters.fa:2:40:15:8:true SLIDINGWINDOW:4:15 LEADING:15 TRAILING:15 MINLEN:36 &\n";
	}
	else
	{
		$trim.="java -jar /proj/b2012111/private/UserDirectory/nima/glob_old/Software/Trimmomatic-0.30/trimmomatic-0.30.jar SE -threads 2 -phred33 -trimlog $lineArr[0].log  $lineArr[3]/$lineArr[1] $for-trimmed.fq ILLUMINACLIP:adapters.fa:2:40:15:8:true SLIDINGWINDOW:4:15 LEADING:15 TRAILING:15 MINLEN:36 &\n";
	}
}
Sbatch_run($trim,"trim",$proj,$cntr*120);
$trim="";
close inSample;
##Recheck quality by FastQC 

my $QC2="module load FastQC\n";
$QC2.="cd $cWD/fastQC-after-trimming/\n";
open(inSample,"Samples.list");
while(<inSample>)
{
	chomp($_);
	@lineArr=split("\t",$_);
	$lineArr[1]=~ m/(.*)\.fastq.gz/;
	$for=$1;
	if($data_type eq "PE")
	{
		$lineArr[2]=~ m/(.*)\.fastq.gz/;
		$rev=$1;
		$QC2.="fastqc --extract --contaminants ~/glob/Contamination/adapters-tab-txt.list -o $lineArr[0]  -t 2 $cWD/Trimmed-reads/$lineArr[0]/$for-paired.fq $cWD/Trimmed-reads/$lineArr[0]/$rev-paired.fq &\n";
	}
	else
	{
		$QC2.="fastqc --extract --contaminants ~/glob/Contamination/adapters-tab-txt.list -o $lineArr[0]  -t 2 $cWD/Trimmed-reads/$lineArr[0]/$for-trimmed.fq &\n";
	}
}
Sbatch_run($QC2,"QC-2",$proj,$cntr*40);
$QC2="";
close inSample;
##Check the quality reports
#open(inSample,"Samples.list");
#while(<inSample>)
#{
#	chomp($_);
#	@lineArr=split("\t",$_);
#		
#}
#close inSample;
##GSNAP
open(inSample,"Samples.list");
while(<inSample>)
{
	chomp($_);
	@lineArr=split("\t",$_);
	$lineArr[1]=~ m/(.*)\.fastq.gz/;
	$for=$1;
	$path_to_ref="";$reference_name="";$reference_fai="";
	@referenceArr=split("\/",$reference);
	for(my $refindex=0;$refindex<scalar(@referenceArr)-1;$refindex++)
	{
		$path_to_ref.="$referenceArr[$refindex]/";
	}
	$reference_name=$referenceArr[-1];
	$reference_fai="$path_to_ref/$referenceArr[-1].fai";
	$reference_fa="$path_to_ref/$referenceArr[-1].fa";
#	$reference=~ m/(.*)\/(.*fa*)/;
#	$path_to_ref=$1;
#	$tmp_reference_name=$2;
#	$tmp_reference_name=~ m/(.*)\.fa*/;
#	$reference_name=$1;
#	print "$path_to_ref/$reference_name.fai";<STDIN>;
	$gsnap.="cd $cWD/GSNAP/$lineArr[0]\nmodule load samtools\n";
	if($vcf ne "" && $data_type eq "PE")
	{
        	$lineArr[2]=~ m/(.*)\.fastq.gz/;
	        $rev=$1;
		$vcf=~ m/(.*)\/(.*).vcf/;
		$path_to_vcf=$1;
		$vcf_name=$2;
		$gsnap.="/proj/b2010051/private/UserDirectories/nima/Tools/gmap-2014-12-23/bin/gsnap -m 0.1 --use-sarray=0 -D $path_to_ref/ -d $reference_name -V $path_to_vcf/ -v $vcf_name -o FR -B 5 -N 1 -n 30 -E 4 --maxsearch=$readLength --nthreads=16 --gmap-mode=pairsearch,terminal,improve -A sam -J 33 -O --quiet-if-excessive --read-group-id=$lineArr[0] --read-group-name=$lineArr[0] --read-group-library=200PE --read-group-platform=Illumina $cWD/Trimmed-reads/$lineArr[0]/$for-paired.fq $cWD/Trimmed-reads/$lineArr[0]/$rev-paired.fq | samtools view -F 4 -q 20 -bSt $reference_fai - >$lineArr[0]-gsnap-snpT.bam \n";
	}
	elsif($vcf eq "" && $data_type eq "PE")
	{
		$lineArr[2]=~ m/(.*)\.fastq.gz/;
                $rev=$1;
		$gsnap.="/proj/b2010051/private/UserDirectories/nima/Tools/gmap-2014-12-23/bin/gsnap -m 0.1 -D $path_to_ref/ -d $reference_name -o FR -B 5 -N 1 -n 30 -E 4 --maxsearch=$readLength --nthreads=16 --gmap-mode=pairsearch,terminal,improve -A sam -J 33 -O --quiet-if-excessive --read-group-id=$lineArr[0] --read-group-name=$lineArr[0] --read-group-library=200PE --read-group-platform=Illumina $cWD/Trimmed-reads/$lineArr[0]/$for-paired.fq $cWD/Trimmed-reads/$lineArr[0]/$rev-paired.fq | samtools view -F 4 -q 20 -bSt $reference_fai - >$lineArr[0]-gsnap.bam \n";
	}
	if($vcf ne "" && $data_type eq "SE")
	{
		$vcf=~ m/(.*)\/(.*).vcf/;
		$path_to_vcf=$1;
		$vcf_name=$2;
		$gsnap.="/proj/b2010051/private/UserDirectories/nima/Tools/gmap-2014-12-23/bin/gsnap -m 0.1 --use-sarray=0 -D $path_to_ref/ -d $reference_name -V $path_to_vcf/ -v $vcf_name -o FR -B 5 -N 1 -n 30 -E 4 --maxsearch=$readLength --nthreads=16 --gmap-mode=pairsearch,terminal,improve -A sam -J 33 -O --quiet-if-excessive --read-group-id=$lineArr[0] --read-group-name=$lineArr[0] --read-group-library=200PE --read-group-platform=Illumina $cWD/Trimmed-reads/$lineArr[0]/$for-trimmed.fq | samtools view -F 4 -q 20 -bSt $reference_fai - >$lineArr[0]-gsnap-snpT.bam \n";
	}
	elsif($vcf eq "" && $data_type eq "SE")
	{
		$lineArr[2]=~ m/(.*)\.fastq.gz/;
                $rev=$1;
		$gsnap.="/proj/b2010051/private/UserDirectories/nima/Tools/gmap-2014-12-23/bin/gsnap -m 0.1 -D $path_to_ref/ -d $reference_name -o FR -B 5 -N 1 -n 30 -E 4 --maxsearch=$readLength --nthreads=16 --gmap-mode=pairsearch,terminal,improve -A sam -J 33 -O --quiet-if-excessive --read-group-id=$lineArr[0] --read-group-name=$lineArr[0] --read-group-library=200PE --read-group-platform=Illumina $cWD/Trimmed-reads/$lineArr[0]/$for-trimmed.fq | samtools view -F 4 -q 20 -bSt $reference_fai - >$lineArr[0]-gsnap.bam \n";
	}
	Sbatch_run($gsnap,"$lineArr[0]-aln",$proj,45*60);
	$gsnap="";
}
close inSample;
##Prepare bam files
open(inSample,"Samples.list");
while(<inSample>)
{
	chomp($_);
	@lineArr=split("\t",$_);
	$prep="module load GATK
module load picard
cd $cWD/GSNAP/$lineArr[0]\n";
	$lineArr[1]=~ m/(.*)\.fastq.gz/;
	$for=$1;
	$lineArr[2]=~ m/(.*)\.fastq.gz/;
	$rev=$1;
	if($vcf ne "")
	{
		$prep.="~/private/Pipelines/SNP-calling/GATK/preparing-bam-files-RNA-seq-V3-removind-unmapped-reads.pl -reference $reference_fa -input $lineArr[0]-gsnap-snpT.bam -sort T -markDup F -recal T -vcf $vcf -recalstat T -rg F -nct 2 -dict F -rm-unmapped F -splitN T >prepare-$lineArr[0].sh
bash prepare-$lineArr[0].sh";
	}
	else
	{
	$prep.="~/private/Pipelines/SNP-calling/GATK/preparing-bam-files-RNA-seq-V3-removind-unmapped-reads.pl -reference $reference_fa -input $lineArr[0]-gsnap-snpT.bam -sort T -markDup F -recal F -recalstat F -rg F -nct 2 -dict F -rm-unmapped F -splitN T >prepare-$lineArr[0].sh
bash prepare-$lineArr[0].sh";
	}
#	print $prep;<STDIN>;
#	Sbatch_run($prep,"$lineArr[0]-prep",$proj,24*60);
	Sbatch_run_core($prep,"$lineArr[0]-prep",$proj,24*60);
	$prep="";
}
close inSample;
##Run QoRTs
open(inSample,"Samples.list");
while(<inSample>)
{
        chomp($_);
        @lineArr=split("\t",$_);
	$lineArr[1]=~ m/(.*)\.fastq.gz/;
	$for=$1;
	$gtf=~ m/(.*)\/(.*)\.gtf/;
	$QC_name=$2;
	$cmd_QoRTs.="cd $cWD/GSNAP/$lineArr[0]
mkdir $lineArr[0]_QC_$QC_name
java -Xmx10G -jar /proj/b2010051/private/UserDirectory/nima/nima/Tools/QoRTs_0.3.18/QoRTs.jar QC --chromSizes $reference_fai --title $lineArr[0] --minMAPQ 20 --stranded --rawfastq $cWD/Trimmed-reads/$lineArr[0]/$for-paired.fq --generatePlots $lineArr[0]-gsnap-snpT.Reordered.sort.bam $gtf $cWD/GSNAP/$lineArr[0]/$lineArr[0]_QC_$QC_name/ ";
	Sbatch_run_core($cmd_QoRTs,"$lineArr[0]-QC",$proj,5*60);
	$cmd_QoRTs="";
}
close inSample;
##Run cufflinks
open(inSample,"Samples.list");
while(<inSample>)
{
	chomp($_);
	@lineArr=split("\t",$_);
	$cufflinks_CMD="module load cufflinks/2.2.1
cd $cWD/cufflinks";
	if($gtf ne "")
	{
		$cufflinks_CMD="cufflinks -p 16 -o $lineArr[0] -u --library-type fr-firststrand -g $gtf $cWD/GSNAP/$lineArr[0]/$lineArr[0]-gsnap-snpT.Reordered.sort.bam";
	}
	else
	{
		$cufflinks_CMD.="cufflinks -p 16 -o $lineArr[0] -b -u --library-type fr-firststrand $cWD/GSNAP/$lineArr[0]/$lineArr[0]-gsnap-snpT.Reordered.sort.bam";
	}
#/proj/b2012206/private/nobackup/annotation/Broad_annotation_2013_10_15/Updated_annotations_dec13/rabbit_annotation/Broad-131213-All-annotations.gtf /proj/b2012206/private/nobackup/UserDirectory/nima/nima/Domestication/GSNAP/$a/$a-gsnap-snpT.Reordered.sort.bam";
	Sbatch_run($cufflinks_CMD,"$lineArr[0]-cufflinks",$proj,15*60);
}
close inSample;
##cuffmerge
#/sw/apps/bioinfo/cufflinks/2.2.1/milou/cuffmerge -o Parents-Ensembl73-cuffmerge-141126 -g /home/nimar/b2012206/private/nobackup/annotation/Rabbit-Ensembl-Genes.gtf -p 16 assembly_list.txt
$cuffmerge_CMD="module load cufflinks/2.2.1
cd $cWD/cufflinks
#mkdir -p cuffmerge-$d
cuffmerge -o cuffmerge-$d -g $gtf -p 16 assembly_list.txt\n";
open(outAssemblyList,">$cWD/GSNAP/assembly_list.txt");
open(inSample,"Samples.list");
while(<inSample>)
{
	chomp($_);
	@lineArr=split("\t",$_);
	$assembly_Files.="$cWD/GSNAP/$lineArr[0]/transcripts.gtf\n";
}
close inSample;
print outAssemblyList $assembly_Files;
Sbatch_run($cuffmerge_CMD,"Sbatch-cuffmerge",$proj,15*60);
$assembly_Files="";


###############Subroutines
sub no_Print_CMD()
{
#	print "CMD: $_[0]\n";
	system ("$_[0]");
}
sub Print_CMD()
{
	print "CMD: $_[0]\n";
	system ("$_[0]");
}
sub Sbatch_run()
{
	$commands=$_[0];
	$step=$_[1];
	$project=$_[2];
	$time=sprintf("%.0f",($_[3]/60));
	open(outSbatch,">$step.script");
	print outSbatch "#!/bin/bash -l
#SBATCH -A $project
#SBATCH -p node
#SBATCH -t $time:00:00
#SBATCH -J $step
#SBATCH --mail-type=all
#SBATCH --mail-user=nimarafati\@gmail.com
$commands
";
	$commands=~ m/.*(\&$)/;
	if($1 eq "&")
	{
		print outSbatch "wait";
	}
}
sub Sbatch_run_core()
{
	$commands=$_[0];
	$step=$_[1];
	$project=$_[2];
	$time=sprintf("%.0f",($_[3]/60));
	open(outSbatch,">$step.script");
	print outSbatch "#!/bin/bash -l
#SBATCH -A $project
#SBATCH -p core -n 2
#SBATCH -t $time:00:00
#SBATCH -J $step
#SBATCH --mail-type=all
#SBATCH --mail-user=nimarafati\@gmail.com
$commands
";
	$commands=~ m/.*(\&$)/;
	if($1 eq "&")
	{
		print outSbatch "wait";
	}
}
