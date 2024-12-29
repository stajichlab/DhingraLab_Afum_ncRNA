#!/usr/bin/bash
#SBATCH -p short --mem 32gb -N 1 -n 16 --out logs/STAR.%a.log -J STAR

FWDEXT=1
REVEXT=2
FASTQEXT=fq.gz

module load star
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi

SAMPLEFILE=samples.csv
INDIR=input
for STRAIN in Af293 A1163
do
    OUTDIR=results/STAR_$STRAIN
    IDX=db/STAR_$STRAIN
    SPECIES=Afumigatus$STRAIN
    GENOME=db/${SPECIES}_Genome.fasta
    GFF=db/${SPECIES}_Genes.gff3
    GTF=db/${SPECIES}_Genes.gtf
    mkdir -p $OUTDIR
    if [ ! -f $GTF ]; then
	grep -P "\texon\t" $GFF | perl -p -e 's/ID=[^;]+;Parent=([^;]+);/gene_id $1/' > $GTF
    fi
    if [ ! -d $IDX ]; then
	STAR --runThreadN $CPU --runMode genomeGenerate --genomeDir $IDX --genomeFastaFiles $GENOME \
	     --sjdbGTFfile $GTF --genomeSAindexNbases 11
    fi
    
    IFS=,
    tail -n +2 $SAMPLEFILE | sed -n ${N}p | while read SAMPLE CONDITION REP READBASE
    do
	OUTNAME=$SAMPLE
	echo "$INDIR/${READBASE}${FWDEXT}.$FASTQEXT $INDIR/${READBASE}_${REVEXT}.$FASTQEXT"
	
	STAR --outSAMstrandField intronMotif --runThreadN $CPU --outMultimapperOrder Random --twopassMode Basic \
	     --genomeDir $IDX --outFileNamePrefix $OUTDIR/$OUTNAME. --readFilesCommand zcat \
	     --readFilesIn $INDIR/${READBASE}${FWDEXT}.$FASTQEXT $INDIR/${READBASE}${REVEXT}.$FASTQEXT
    done
    unset IFS
done
