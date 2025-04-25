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
    GO=${SPECIES}_GO.gaf
    GOCURATED=${SPECIES}_Curated_GO.gaf

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
        curl -o $DB/$GO.gz $URL/$SPECIES/gaf/FungiDB-${VERSION}_${GO}.gz
	gunzip $DB/$GO.gz
    fi
    
    if [ ! -f $DB/${GOCURATED} ]; then
        curl -o $DB/${GOCURATED}.gz $URL/$SPECIES/gaf/FungiDB-${VERSION}_${GOCURATED}.gz
	gunzip $DB/$GOCURATED.gz
    fi
done

if [ ! -f $DB/goslim_aspergillus.obo ]; then
	curl -o $DB/goslim_aspergillus.obo https://current.geneontology.org/ontology/subsets/goslim_aspergillus.obo
fi
if [ ! -f $DB/go.obo ]; then
	curl -L -o $DB/go.obo https://current.geneontology.org/ontology/go.obo
fi
