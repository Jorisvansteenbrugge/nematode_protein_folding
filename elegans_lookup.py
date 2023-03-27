#!/usr/bin/env python3

import pandas as pd
from jorispython import get_uniprot_gene_description

elegans_comparisons = pd.read_csv("celegans/All_elegans_pallida_bonferroni.csv", sep = ';')

elegans_comparisons['Query'] = [ x.split('/')[-1].replace(".pdb","") for x in elegans_comparisons['Query']]
elegans_comparisons['Subject'] = [x.split('/')[-1].replace(".pdb","") for x in elegans_comparisons['Subject']]

#pvalue commas to periods
elegans_comparisons["P.adjusted"]= pd.to_numeric([x.replace(",",".") for x in elegans_comparisons["P.adjusted"]])

# uit unknown genome.py
gpal_unknown_genes = pd.read_csv("knowable_genome/gpal_unannotated_gene_ids.csv", sep = ';')["Unannotated_gene"]

gpal_unknown_genes = [x.replace(".t1","") for x in gpal_unknown_genes ]

elegans_comparisons_significant = elegans_comparisons[elegans_comparisons["P.adjusted"] <= 0.05]

gp_unknown_genes_match = elegans_comparisons_significant.query("Subject==@gpal_unknown_genes")
ids_clean = list([ele_id.split("-")[1] for ele_id in gp_unknown_genes_match["Query"]])

gp_unknown_genes_match_clean = gp_unknown_genes_match
gp_unknown_genes_match_clean["Query"] = ids_clean

gp_unknown_genes_match_clean.to_csv("celegans/pallida_elegans_hits.csv")

# create network file
elegans_id_annotations = pd.read_csv("celegans/elegans_functional_annotations.tsv", sep ='\t')

with open("celegans/pallida_elegans_hits_annotated.csv",'w') as network_file:
    network_file.write("Gpal_ID\tUniprot_ID\tProtein Name\t Gene Name\t Wormbase ref\tSpec_A\tSpec_B\n")
    for gpal_id in gpal_unknown_genes:
        for elegans_match in elegans_comparisons_significant.query("Subject==@gpal_id")["Query"]:
            elegans_match_id = elegans_match.split('-')[1]
            elegans_anno = elegans_id_annotations[
                (elegans_id_annotations['Entry'] ==elegans_match_id)][ [
                    "From", "Protein names", "Gene Names", "WormBase"] ]
            if len(elegans_anno.index) == 0:
                annotation = "NA"
            else:
                annotation = "\t".join([str(x) for x in elegans_anno.iloc[0]])
            network_file.write(f"{gpal_id}\t{annotation}\tGpal\tCele\n")
