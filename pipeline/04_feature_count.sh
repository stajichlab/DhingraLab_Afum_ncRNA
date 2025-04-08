#!/usr/bin/bash -l
#SBATCH -p epyc  -N 1 -n 1 -c 24 --mem 64gb --out logs/subread.log --time 4:00:00

module load subread
CPU=24

featureCounts -g gene_id -G db/AfumigatusA1163_Genome.fasta -s 0 \
	-a db/AfumigatusA1163_Genes.gtf --tmpDir $SCRATCH -o results/read_count_A1163_sorted.tsv -F GTF \
	-T $CPU -M results/STAR_A1163_sort/*.bam -p

featureCounts -g gene_id -G db/AfumigatusAf293_Genome.fasta -s 0 \
	-a db/AfumigatusAf293_Genes.gtf --tmpDir $SCRATCH -o results/read_count_Af293_sorted.tsv -F GTF \
	-T $CPU -M results/STAR_Af293_sort/*.bam -p
