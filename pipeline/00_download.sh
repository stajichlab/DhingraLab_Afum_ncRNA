#!/usr/bin/bash
#SBATCH -p short

# FIX THIS TO SPECIES YOU WANT FROM FUNGIDB 
# OR REWRITE THIS SCRIPT TO DOWNLOAD WHAT YOU WANT
VERSION=68
DB=db
mkdir -p $DB

for SPECIES in AfumigatusAf293 AfumigatusA1163
do
    URL=https://fungidb.org/common/downloads/release-${VERSION}
    MRNA=${SPECIES}_mRNA.fasta
    GENOME=${SPECIES}_Genome.fasta
    GFF=${SPECIES}_Genes.gff3
    GO=${SPECIES}_GO.gaf.gz
    GOCURATED=${SPECIES}_Curated_GO.gaf.gz

    if [ ! -f $DB/$MRNA ]; then
        curl -o $DB/$MRNA $URL/$SPECIES/fasta/data/FungiDB-${VERSION}_${SPECIES}_AnnotatedTranscripts.fasta 
    fi
    
    if [ ! -f $DB/$GENOME ]; then
        curl -o $DB/$GENOME $URL/$SPECIES/fasta/data/FungiDB-${VERSION}_${SPECIES}_Genome.fasta
    fi
    
    if [ ! -f $DB/$GFF ]; then
        curl -o $DB/$GFF $URL/$SPECIES/gff/data/FungiDB-${VERSION}_${SPECIES}.gff
    fi
    if [ ! -f $DB/${GO} ]; then
        curl -o $DB/$GO $URL/$SPECIES/gaf/data/FungiDB-${VERSION}_${GO}
    fi

    if [ ! -f $DB/${GOCURATED} ]; then
        curl -o $DB/${GOCURATED} $URL/$SPECIES/gaf/data/FungiDB-${VERSION}_${GOCURATED}
    fi


done
