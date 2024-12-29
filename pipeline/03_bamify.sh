#!/usr/bin/bash -l
#SBATCH -p epyc -c 24 -N 1 -n 1 --mem 64gb --out logs/make_bam.log

CPU=4
module load samtools
FOLDER=results
for STRAIN in A1163 Af293
do
	parallel -j 6 samtools view --threads $CPU -O BAM -o {.}.bam {} ::: $(find $FOLDER/STAR_$STRAIN -name "*.sam")
done
