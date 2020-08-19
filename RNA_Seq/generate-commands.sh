#!/usr/bash
###################################################################
#0. Resources                                                     #
###################################################################
PROJ_ID='snic2020-15-15'
PROJECT_DIR='/crex/proj/snic2020-16-21/nobackup/private/SMS_5221_20_SounDrivenBiotechnology/'
CODE=$PROJECT_DIR'code/'
SAMPLES=$CODE'Samples_reads.txt'
INTERMEDIATE_DIR=$PROJECT_DIR'intermediate/'
RESULTS_DIR=$PROJECT_DIR'results/'
RAW_READS=$PROJECT_DIR'data/raw_data/'
TRIMMED_READS=$PROJECT_DIR'data/raw_data/Trimmed-reads/'
STAR=$INTERMEDIATE_DIR'STAR_chromosome_level_Trimmed/'
STAR_QC=$INTERMEDIATE_DIR'STAR_chromosome_level_Trimmed_QC/'
GENOME_DIR=$PROJECT_DIR'data/meta_data/reference/STARIndex/'
GFF_FILE=$PROJECT_DIR'data/meta_data/annotation/annotation.gff'
GTF_FILE=$PROJECT_DIR'data/meta_data/annotation//annotation.gtf_new'
GENOME_SIZE=$PROJECT_DIR'data/meta_data/reference/genome.fai'
FEATURECOUNTS_DIR=$PROJECT_DIR'intermediate/featureCounts/STAR_chromosome_level_Trimmed/'

mkdir $TRIMMED_READS $STAR $STAR_QC $FEATURECOUNTS_DIR
THREADS=20

SBATCH="#!/bin/bash -l
#SBATCH -A $PROJ_ID
#SBATCH --mail-type=all
#SBATCH --mail-user=nimarafati@gmail.com"

###################################################################
# FastQC                                                     	  #
###################################################################
rm -rf run_fastqc.txt
while read -r sample R1 R2
do
	echo "$SBATCH
#SBATCH -J ${sample}_FastQC
#SBATCH -p core -n 4
#SBATCH -t 2:00:00
module load bioinfo-tools FastQC

cd \$SNIC_TMP
cp $RAW_READS/$sample/*fastq.gz  \$SNIC_TMP/ 
mkdir  ${sample}_QC
mkdir -p $INTERMEDIATE_DIR/FastQC/
fastqc -t 4 -o ${sample}_QC -f fastq -t 4 *fastq.gz
mv ${sample}_QC $INTERMEDIATE_DIR/FastQC/" >Sbatch_fastqc_${sample}.script
echo "sbatch Sbatch_fastqc_${sample}.script" >>run_fastqc.txt
done<$SAMPLES
#Samples_reads.txt

###################################################################
# Trimming                                                        #
###################################################################
rm -f run_trim.log
THREADS=2
while read -r sample R1 R2
do
	echo "$SBATCH
#SBATCH -J ${sample}_trim
#SBATCH -p core -n $THREADS
#SBATCH -t 6:00:00
module load bioinfo-tools trimmomatic
mkdir $TRIMMED_READS/${sample}

cp $RAW_READS/$sample/*fastq.gz \$SNIC_TMP/
cd \$SNIC_TMP
cp /home/nimar/glob_old/Contamination/adapters.fa adapters.fa 
trimmomatic PE -threads $THREADS -phred33 -trimlog ${sample}.log \
$R1 $R2 \
${sample}_R1_paired.fq ${sample}_R1_unpaired.fq \
${sample}_R2_paired.fq ${sample}_R2_unpaired.fq \
ILLUMINACLIP:adapters.fa:2:40:15:8:true SLIDINGWINDOW:4:15 LEADING:10 TRAILING:10 MINLEN:36 HEADCROP:10
#Compress paired reads
gzip ${sample}_R1_paired.fq &
gzip ${sample}_R2_paired.fq &
wait 
#Transfer only paired reads
mv ${sample}_R1_paired.fq.gz $TRIMMED_READS/${sample} &
mv ${sample}_R2_paired.fq.gz $TRIMMED_READS/${sample} &
wait" >Sbatch_trim_${sample}.script
	echo "sbatch Sbatch_trim_${sample}.script" >>run_trim.log
done<$SAMPLES

###################################################################
# Align the reads to genome and Index bam files                   #
###################################################################
rm -f run_align.sh
while read -r sample R1 R2
do
        echo "$SBATCH
#SBATCH -p node
#SBATCH -t 1:00:00
#SBATCH -J ${sample}_align

module load bioinfo-tools star samtools
cd \$SNIC_TMP
cp $TRIMMED_READS/*$sample*/$R1 \$SNIC_TMP &
cp $TRIMMED_READS/*$sample*/$R2 \$SNIC_TMP &
wait

#Alignign the reads
mkdir -p $STAR/$sample
STAR --genomeDir $GENOME_DIR \
--sjdbGTFfile $GTF_FILE \
--readFilesIn $R1 $R2 \
--runThreadN  $THREADS \
--twopassMode Basic \
--outWigType bedGraph \
--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 64324509440 \
--readFilesCommand zcat \
--runDirPerm All_RWX  \
--quantMode GeneCounts \
--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 20 \
--outFileNamePrefix $sample --outSAMattrRGline ID:$sample 'SM:${sample}' 

#Simplify name of the bam and count files 
mv ${sample}Aligned.sortedByCoord.out.bam ${sample}.sort.bam

#Indexing bam files
samtools index -@ $THREADS ${sample}.sort.bam
ln -s  ${sample}.sort.bam.bai ${sample}.sort.bai
mv ${sample}* $STAR/$sample/
cd $CODE
sbatch Sbatch_featurecounts_${sample}.script" >Sbatch_align_${sample}.script
echo "sbatch Sbatch_align_${sample}.script" >>run_align.sh
done<$SAMPLES
###################################################################
# QC by QoRTs, qualimap                                           #
###################################################################
#Runnign QoRTs
rm -f run_QoRTs.sh
THREADS=20
while read -r sample R1 R2
do
        echo "$SBATCH
#SBATCH -p node
#SBATCH -t 1:30:00
#SBATCH -J QoRTS_${sample}
#module load bioinfo-tools QualiMap

java -Xmx50G -jar ~/git/QoRTs/QoRTs.jar QC --prefilterImproperPairs --generatePlots --stranded --generatePlots --numThreads $THREADS --maxReadLength $READ_LENGTH  --minMAPQ 20 --title ${sample} \
--addFunctions JunctionCalcs,writeJunctionSeqCounts,writeSpliceExon,writeNovelSplices,annotatedSpliceExonCounts,makeAllBrowserTracks,FPKM,writeGeneCounts,writeDESeq,writeDEXSeq,writeGeneBody,writeGenewiseGeneBody,calcDetailedGeneCounts,makeJunctionBed,writeDocs,writeGeneBodyIv --chromSizes $GENOME_SIZE $STAR/$sample/${sample}.sort.bam \
$GTF_FILE $STAR_QC${sample}_QoRTs/ 
#qualimap bamqc -nt 20 -bam $STAR/$sample/${sample}.sort.bam -c -gff $GFF_FILE -outdir $STAR_QC${sample}_qualimap -p strand-specific-reverse  -sd " >Sbatch_QoRTs_${sample}.script
echo "sbatch Sbatch_QoRTs_${sample}.script" >>run_QoRTs.sh
done<$SAMPLES


###################################################################
# featurecounts                                                   #
###################################################################
rm -f run_featureCounts.sh
while read -r sample R1 R2
do
        echo "$SBATCH
#SBATCH -p node
#SBATCH -t 00:30:00
#SBATCH -J ${sample}_FC
cd \$SNIC_TMP
cp $STAR/$sample/${sample}.sort.bam \$SNIC_TMP
featureCounts_path=\"$FEATURECOUNTS_DIR/${sample}\"
output=\"count-s-2\"
annotation=\"$GTF_FILE\"
mkdir -p \$featureCounts_path
~/glob/Software/subread-2.0.0-source/bin/featureCounts -s 2 -t exon -g locus_tag -F GTF -C -T $THREADS -p -B --primary -Q 20 \
-o \$output \
-a \$annotation \
${sample}.sort.bam 
mv \$output* \$featureCounts_path" >Sbatch_featurecounts_${sample}.script
echo "sbatch Sbatch_featurecounts_${sample}.script " >>run_featureCounts.sh
done<$SAMPLES
