#!/usr/bin/env python3

"""Get the knowable genome - i.e. functionally annotate genes based on other genes that have a
similar shape
"""


import pandas as pd
import logging
import venn # local py file

logging.basicConfig(level=logging.INFO)

gp_unkown = []
gp_knowable = []

mc_unkown = []
mc_knowable = []

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
    if spprefix == "Gpal":
        knowable = gp_knowable
        unkown = gp_unkown
    elif spprefix == "Mchit":
        knowable = mc_knowable
        unkown = mc_unkown

    if False in annotated:
        if True in annotated:
            knowable.append(str(cluster))
        else:
            unkown.append(str(cluster))


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


labels = venn.get_labels([gp_knowable,
                          gp_unkown,
                          mc_knowable,
                          mc_unkown])
fig, ax = venn.venn4(labels, names = ["Gp Knowable",
                                      "Gp Unknown",
                                      "Mc Knowable",
                                      "Mc Unknown"])
fig.show()
