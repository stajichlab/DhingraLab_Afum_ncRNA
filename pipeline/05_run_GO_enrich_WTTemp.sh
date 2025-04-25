#!/usr/bin/bash -l
#SBATCH -p short -c 2 --mem 2gb --job-name GO_enrich --out logs/GO_enrich_WTvDel.out

for file in $(ls reports/pairwise_exp_direction/venncount_*_names.tsv)
do
    bash scripts/GoEnrichmentOfGeneIDs.sh $file
    if [ $? -ne 0 ]; then
        echo "Error running GoEnrichmentOfGeneIDs.sh on $file"
        exit 1
    fi
    echo "GO enrichment analysis completed for $file"
done
