#!/usr/bin/bash -l
#SBATCH -p short -c 16 --mem 24gb --out logs/ncRNA_expression.log

module load subread

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

mkdir -p results

featureCounts -G ../db/AfumigatusAf293_Genome.fasta -s 0 \
	-a ncRNA.saf --tmpDir $SCRATCH -o results/WT_ncRNA_read_count_Af293 -F SAF \
	-T $CPU -M ../results/STAR_Af293/WT*.bam -p 


featureCounts -G ../db/AfumigatusAf293_Genome.fasta -s 0 \
	-a ncRNA.saf --tmpDir $SCRATCH -o results/DEL_ncRNA_read_count_Af293 -F SAF \
	-T $CPU -M ../results/STAR_Af293/DEL*.bam -p
