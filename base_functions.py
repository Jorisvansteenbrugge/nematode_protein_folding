#!/usr/bin/env python3
import pandas ad pd

cluster_df = pd.read_csv("cytoscape_clusters_new.csv")
edges = pd.read_csv("cytoscape_edges.csv", sep =';')


def get_gene_descriptions(list: gene_s):
    return cluster_df.query('id in @gene_s')
