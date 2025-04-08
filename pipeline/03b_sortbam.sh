#!/usr/bin/bash -l
#SBATCH -c 24 -N 1 -n 1 --mem 64gb --out logs/make_bam_sort.log

CPU=24
module load samtools
FOLDER=results
	#Af293
for STRAIN in Af293 A1163
do
	mkdir -p $FOLDER/STAR_${STRAIN}_sort
	for file in $(find $FOLDER/STAR_$STRAIN -name "*.bam")
	do
		OUT=$(basename $file .Aligned.out.bam)
		if [ ! -f $FOLDER/STAR_${STRAIN}_sort/$OUT.bam.csi ]; then
			samtools sort -T $SCRATCH/$OUT --threads $CPU --write-index -o $FOLDER/STAR_${STRAIN}_sort/$OUT.bam $file
		fi
	done
done

