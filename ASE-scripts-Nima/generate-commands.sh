#!/usr/bash
while read -r a b c;
do
	#keep the best reads from alignment on each parental genome (../Samples.list)
#	echo "./keep_best.pl $a /proj/b2011054/private/Users/nima/tmpRabbit/oryCun2_assembly_files/oryCun2-reference-GSNAP/oryCun2-GSNAP.fa &"
	echo "samtools view -bS -t /proj/b2011054/private/Users/nima/tmpRabbit/oryCun2_assembly_files/oryCun2-reference-GSNAP/oryCun2-GSNAP.fa.fai $a.sam >$a.bam &" 
done<sam.list
#../Samples.list
