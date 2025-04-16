#!/usr/bin/env python3

import csv
import os

infolder = "reports"
outfolder = "reports/venngroups"



def get_pairwise(infile):
    with open(infile, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        genes = {}
        for row in reader:
            # Check if the row contains the expected number of columns
            if len(row) != len(header):
                print(f"Skipping row with unexpected number of columns: {row}")
                continue
            gene = row[0]
            foldchange = row[2]
            if gene not in genes:
                genes[gene] = foldchange
            else:
                print(f"Duplicate gene found: {gene}. Skipping this gene.")
                continue
        return genes

def main(inf=infolder,outf=outfolder):
    if not os.path.exists(outf):
        os.makedirs(outf)
        print(f"Created folder: {outfolder}")
    groupset = {}
    for reportfile in os.listdir(inf):
        if reportfile.endswith(".csv") and "_vs_" in reportfile:
            infile = os.path.join(inf, reportfile)
            (groupname,compare1,ignore,compare2) = reportfile.split('.csv')[0].split("_")
            if groupname == "all":
                continue
            compare = compare1 + "_v_" + compare2
            print(f"Processing file: {reportfile}, group: {groupname}, comparison: {compare}")

            groups = get_pairwise(infile)
            if groupname not in groupset:
                groupset[groupname] = {}
            for gene in groups:
                if gene not in groupset[groupname]:
                    groupset[groupname][gene] = {}
                groupset[groupname][gene][compare] = float(groups[gene])   # storing fold change
            # Create the output file name    
    for groupname, genes in groupset.items():            
            # Write the groups to the output file
            gene_classification = {
                                'up': {}, 
                                'down': {},
                                'all': {}
                                }
            for gene in genes:
                classification = {'up': [], 'down': [], 'all': []}
                
                for compare in genes[gene]:
                    #print(f"Processing gene: {gene}, comparison: {compare}, fold change: {genes[gene][compare]}")
                    classification['all'].append(compare)
                    if genes[gene][compare] > 0:
                        classification['up'].append(compare)
                    else:
                        classification['down'].append(compare)
                                
                for direction in classification:
                    if len(classification[direction]) > 0:
                        class_string = ",".join(sorted(classification[direction]))
                        if class_string not in gene_classification[direction]:
                            gene_classification[direction][class_string] = set()
                        gene_classification[direction][class_string].add(gene)

            # Create the output file name for all groups
            for direction in ['up', 'down', 'all']:
                outfilect = os.path.join(outf, f"venncount_{groupname}_{direction}_count.tsv")
                outfilenames = os.path.join(outf, f"venncount_{groupname}_{direction}_names.tsv")
                with open(outfilect, 'w', newline='') as ctfh, open(outfilenames, 'w', newline='') as namefh:
                    writer = csv.writer(namefh, delimiter="\t")
                    writer.writerow(["Gene", "Classification"])
                    countwriter = csv.writer(ctfh, delimiter="\t")
                    countwriter.writerow(["Classification","Count"])
                    for classtype,class_set in gene_classification[direction].items():
                        countwriter.writerow([classtype,len(class_set)])
                        for gene in class_set:
                            writer.writerow([gene, classtype])

main()