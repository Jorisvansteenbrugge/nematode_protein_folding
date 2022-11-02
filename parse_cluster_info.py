import pandas as pd
import seaborn as sns

print()
df = pd.read_csv("cytoscape_clusters_new.csv")
clusters = df["__fastGreedyCluster"]


def get_identities(edges, nodes):
    subset_edges = edges.loc[edges['source'].isin(nodes) | edges['target'].isin(nodes)]
    return subset_edges[['Identity','P.adjusted']]

def get_cluster_identities():
    output = []

    edges = pd.read_csv("cytoscape_edges.csv", sep=';')

    unique_clusters = list(set(clusters))

    
    for cluster in unique_clusters:
        print(f"Cluster: {cluster}")
        nodes_cluster_c = df.loc[df['__fastGreedyCluster'] == cluster]
        nodes_c = nodes_cluster_c['name']
        
        identities_c = get_identities(edges, list(nodes_c))
        identities_pvalue_rows = identities_c
        identities_pvalue_rows['Cluster']=f"Cl{cluster}"
        output.append(identities_pvalue_rows)

    output = pd.concat(output)

    return output
