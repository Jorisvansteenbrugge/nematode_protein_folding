#!/usr/bin/env python3

import pandas as pd
from jorispython import get_uniprot_gene_description

elegans_comparisons = pd.read_csv("celegans/All_elegans_pallida_bonferroni.csv", sep = ';')

elegans_comparisons['Query'] = [ x.split('/')[-1].replace(".pdb","") for x in elegans_comparisons['Query']]
elegans_comparisons['Subject'] = [x.split('/')[-1].replace(".pdb","") for x in elegans_comparisons['Subject']]

#pvalue comma's to .
elegans_comparisons["P.adjusted"]= pd.to_numeric([x.replace(",",".") for x in elegans_comparisons["P.adjusted"]])



# uit unknown genome.py
gpal_unknown_genes = ("Gpal_D383_g17406.t1","Gpal_D383_g13571.t1","Gpal_D383_g08491.t1",
                     "Gpal_D383_g07960.t1","Gpal_D383_g00459.t1","Gpal_D383_g04885.t1",
                     "Gpal_D383_g01745.t1","Gpal_D383_g14216.t1","Gpal_D383_g06928.t1",
                     "Gpal_D383_g05579.t1","Gpal_D383_g10018.t1","Gpal_D383_g12355.t1",
                     "Gpal_D383_g11164.t1")

gpal_unknown_genes = [x.replace(".t1","") for x in gpal_unknown_genes ]



elegans_comparisons_significant = elegans_comparisons[elegans_comparisons["P.adjusted"] <= 0.05]

gp_unknown_genes_match = elegans_comparisons_significant.query("Subject==@gpal_unknown_genes")
ids_clean = list(set([ele_id.split("-")[1] for ele_id in gp_unknown_genes_match["Query"]]))



# create network file
elegans_id_annotations = pd.read_csv("celegans/elegans_mart_export.csv", sep =',')
with open("elegans_functions_network.csv",'w') as network_file:
    for gpal_id in gpal_unknown_genes:
        for elegans_match in elegans_comparisons_significant.query("Subject==@gpal_id")["Query"]:
            elegans_match_id = elegans_match.split('-')[1]
            elegans_anno = elegans_id_annotations[(elegans_id_annotations['UniProtKB Gene Name ID'] ==elegans_match_id)]
            if len(elegans_anno.index) == 0:
                annotation = "NA"
            else:
                annotation = ",".join([str(x) for x in elegans_anno.iloc[0]])
            network_file.write(f"{gpal_id};{elegans_match};{annotation};Gpal;Cele\n")
