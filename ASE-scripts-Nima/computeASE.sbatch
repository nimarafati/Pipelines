#!/bin/bash -l
#SBATCH -A b2013097
#SBATCH -p core
#SBATCH -t 24:00:00
#SBATCH -J ASE
#SBATCH -o ASE.log
#SBATCH -e ASE.log

prefix="."
gatk="java -Xmx23g -jar /sw/apps/bioinfo/GATK/3.2.2/GenomeAnalysisTK.jar"
ref="/proj/b2013097/private/data/Reference/Galgal4.fa"

#cat pass_hz_exon.vcf | awk -v OFS='\t' '{print $1,$2,$2+1}' > regions.bed

# Use file without IndDPFilter and without repeatmasker regions
#zcat $prefix/alt_ref/HC.noIndDP.vcf.gz > $prefix/offspring/ASE/HC.noIndDP.vcf

#$gatk -T SelectVariants -R $ref -o $prefix/offspring/ASE/HC.filtered.SNP.vcf --variant $prefix/offspring/ASE/HC.noIndDP.vcf --excludeFiltered -selectType SNP

#bedtools intersect -a $prefix/offspring/ASE/HC.filtered.SNP.vcf -b /proj/b2013097/private/data/Reference/Gallus_gallus.Galgal4.73.gtf -u > $prefix/offspring/ASE/HC.filtered.SNP.exon.vcf

#cat $prefix/offspring/ASE/HC.filtered.SNP.exon.vcf | grep -v "#" | grep "0/0" | grep "1/1" > $prefix/offspring/ASE/HC.filtered.SNP.exon.hz.vcf

#ln -s /proj/b2013097/nobackup/private/UserDirectories/susanneB/offspring/ASE/HC.filtered.SNP.exon.hz.vcf variants.vcf

mkdir -p out
i=0

parents=`ls -d 1*`
for parent in $parents
do
	echo Parent: $parent
	mkdir -p out/$parent
	bams=`ls $parent/*.sort.bam`
	for bam in $bams
	do
		echo $bam
	   	sample=`echo $bam | cut -d'/' -f 2 | cut -d'.' -f 1`
	   	date=`echo $bam | cut -d'/' -f 2 | cut -d'.' -f 2`
		if [ -f "out/$parent/$sample.$date.out" ]; then
			echo $bam is already done. Skipping...
		else
			i=$((i+1))
			perl computeASE.pl variants.vcf $bam > out/$parent/$sample.$date.out 2> out/$parent/$sample.$date.err &
			echo $i
			if [ $i -ge 4 ]; then
				wait;
				i=0
			fi
		fi
	done
done
