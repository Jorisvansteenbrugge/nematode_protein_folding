#!/usr/bin/env python3

"""DEPRECATED SCRIPT"""



"""Get the knowable genome - i.e. functionally annotate genes based on other genes that have a
similar shape
"""


import pandas as pd
import logging
import venn # local py file
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.INFO)

gp_unkown   = []
gp_knowable = []
gp_knowable_genes = []
gp_unknown_genes   = []

mc_unkown   = []
mc_knowable = []
mc_knowable_genes = []
mc_unknown_genes  = []

def get_genes_cluster(cluster_name, network_table, species_marker):
    sub = network_table[network_table["__fastGreedyCluster"] == cluster_name]
    return list(sub[sub["species"] == species_marker]["shared name"])

def clean_gene_names(cluster, network_table, species_marker):
    genes = get_genes_cluster(cluster, network_table, species_marker)
    return [gene.replace(".pdb", ".t1") for gene in genes]

def _has_interGO(value):
    return value not in ("no IPS match", "no GO terms")

def _has_description(value):
    if "-NA-" in value:
        return False
    else:
        return True

def is_annotated(row, t):
    description = _has_description(row["Description"])
    interGO = _has_interGO(row["InterPro GO IDs"])
    return True in (description, interGO)

def check_annotations(fun_anno, spprefix, annotation_table):
    annotated = list(
        fun_anno.apply(is_annotated,
                            axis=1,
                            t=gp_gene_annotations))

    unannotated = [not x for x in annotated]
    unannotated_gene_ids = list(fun_anno[unannotated]["SeqName"])

    if spprefix == "Gpal":
        knowable       = gp_knowable
        unkown         = gp_unkown
        knowable_genes = gp_knowable_genes
        unknown_genes  = gp_unknown_genes

    elif spprefix == "Mchit":
        knowable       = mc_knowable
        unkown         = mc_unkown
        knowable_genes = mc_knowable_genes
        unknown_genes  = mc_unknown_genes

    if False in annotated:
        if True in annotated:
            knowable.append(str(cluster))
            knowable_genes += unannotated_gene_ids
        else:
            unkown.append(str(cluster))
            unknown_genes += unannotated_gene_ids


def get_knowable_gpal():
    return gp_knowable

network_table = pd.read_csv("cytoscape_clusters.csv",
                            sep=',')

gp_gene_annotations = pd.read_csv("/mnt/nemahomes/steen176/Genomes/Annotations/G_pallida_annotations.csv",
                                   sep="\t")

mc_gene_annotations = pd.read_csv("/mnt/dataHDD/scratch/chitwoodi_genomes/functional_annotation2.csv",
                                  sep='\t')

clusters = list(set(network_table["__fastGreedyCluster"]))
for cluster in clusters:
    gp_genes = clean_gene_names(cluster, network_table, "Gpal")
    gp_functional = gp_gene_annotations.query("SeqName.isin(@gp_genes)")

    mc_genes = clean_gene_names(cluster, network_table, "Mchit")
    mc_functional = mc_gene_annotations.query("SeqName.isin(@mc_genes)")

    check_annotations(gp_functional, "Gpal", gp_gene_annotations)
    check_annotations(mc_functional, "Mchit", mc_gene_annotations)

logging.info(f"Gp knowable: {len(gp_knowable)}")
logging.info(f"Gp unkown: {len(gp_unkown)}")

logging.info(f"Gp knowable clusters: {', '.join(gp_knowable)}")
logging.info(f"Gp unkown clusters: {', '.join(gp_unkown)}")

logging.info(f"Mc knowable: {len(mc_knowable)}")
logging.info(f"Mc unkown: {len(mc_unkown)}")

logging.info(f"Mc knowable clusters: {', '.join(mc_knowable)}")
logging.info(f"Mc unkown clusters: {', '.join(mc_unkown)}")

print(gp_knowable_genes)

labels = venn.get_labels([gp_knowable, mc_knowable])

fig1,ax1 = venn.venn2(labels, names = ["Gp Knowable", "Mc Knowable"])

fig1.savefig("clusters_knowable_venn.png")
fig1.clf()

labels2 = venn.get_labels([gp_unkown, mc_unkown])
fig2,ax2 = venn.venn2(labels2, names = ["GP Unknown", "Mc Unknown"])

fig2.savefig("clusters_unknown_venn.png")

print("Gpal genes: {}".format(network_table.query("species=='Gpal'").shape[0]))
print("Mchit genes: {}".format(network_table.query("species=='Mchit'").shape[0]))

elegans_comparisons = pd.read_csv("celegans/All_elegans_pallida", sep = '\t')

elegans_comparisons['Query']   = [ x.split('/')[-1] for x in elegans_comparisons['Query'] ]
elegans_comparisons['Subject'] = [ x.split('/')[-1] for x in elegans_comparisons['Subject'] ]

for u_cluster in gp_unkown:
    genes = clean_gene_names(int(u_cluster), network_table, "Gpal")
    break

"""
// Enter Flows between Nodes, like this:
//         Source [AMOUNT] Target

G pallida Proteins [406] Known
G pallida Proteins [296] Knowable
G pallida Proteins [13] Unkown

M chitwoodi Proteins [53] Knowable
M chitwoodi Proteins [3] Unkown
M chitwoodi Proteins [624] Known
"""
