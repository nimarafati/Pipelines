mkdir -p bams_sort
files=`ls [45]*.sam`
for file in $files
do
	name=`echo ${file%.sam}`
	echo $name

   if [ -e "bams_sort/$name.sort.bam" ]; then
      echo NOTE: $file already sorted. Skipping...
   else
      # Avoid sorting files that may be created right now
      if test `find "$file" -mmin +5`; then
			samtools sort -o bams_sort/$name.sort.bam -T /tmp_fast/ $file
			samtools index bams_sort/$name.sort.bam
		else
			echo Too new. Skipping...
		fi
	fi
done
