import pandas as pd
import seaborn as sns
import scipy.stats as stats

from goenrich import *
from collections import Counter

df = pd.read_csv("cytoscape_clusters_new.csv")
edges = pd.read_csv("cytoscape_edges.csv", sep=';')

clusters = df["__fastGreedyCluster"]

GO_ID_background = [extract_go_ids(id_string) for id_string in df["GO.IDs"]]
GO_ID_background = [clean_go_term(go) for go in golist_flattern(GO_ID_background)]



def mean(values):
    values = [float(x.replace(',','.')) for x in values]
    return sum(values)/len(values)

def clean_go_term(go):
    return ":".join(go.split(":")[1:])

def get_cluster_identities():
    output = []

    unique_clusters = list(set(clusters))
    
    for cluster in unique_clusters:
        print(f"Cluster: {cluster}")
        nodes_cluster_c = df.loc[df['__fastGreedyCluster'] == cluster]
        nodes_c = nodes_cluster_c['name']
        
        identities_c = get_identities(edges, list(nodes_c))
        identities_pvalue_rows = identities_c
        identities_pvalue_rows['Cluster'] = f"Cl{cluster}"
        output.append(identities_pvalue_rows)

    output = pd.concat(output)

    return output

def get_enriched(current_terms):
    import obonet
    from goatools.obo_parser import GODag
    from goatools.base import download_go_basic_obo

    # Download the OBO file if it's not already downloaded
    obo_file = 'go-basic.obo'
    if not os.path.exists(obo_file):
        download_go_basic_obo(obo_file)

    # Load the OBO file
    graph = obonet.read_obo(obo_file)

    # Define your lists
    annotation_list = current_terms
    background_list = ["GO:0005737", "GO:0005829", "GO:0006412", "GO:0000002", "GO:0000003", "GO:0000006", "GO:0000012", "GO:0000014", "GO:0000015"]

    # Create a dictionary that maps each GO term ID to its name
    id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}

    # Create a GODag object for the OBO file
    godag = GODag(obo_file)

    # Define the background set
    background_set = set(background_list)

    # Define the study set
    study_set = set(annotation_list)

    # Calculate the number of genes in the study set and background set annotated with each GO term
    go_counts = {}
    for term in background_list:
        annotated_study_genes = len(study_set.intersection(godag.get(term).get_all_children()))
        annotated_background_genes = len(background_set.intersection(godag.get(term).get_all_children()))
        go_counts[term] = {'study': annotated_study_genes, 'background': annotated_background_genes}

    # Perform a Fisher's exact test to determine enriched GO terms
    enriched_terms = []
    for term in go_counts:
        a = go_counts[term]['study']
        b = len(study_set) - a
        c = go_counts[term]['background'] - a
        d = len(background_set) - len(study_set) - c
        oddsratio, pvalue = stats.fisher_exact([[a, b], [c, d]])
        if pvalue < 0.05:
            enriched_terms.append({'id': term, 'name': id_to_name[term], 'pvalue': pvalue})

    # Print the results
    if len(enriched_terms) == 0:
        print('No enriched terms found.')
    else:
        print('Enriched GO terms:')
        for term in enriched_terms:
            print(term['id'], term['name'], term['pvalue'])




def extract_go_ids(id_string):
    """Ids can contain multiple GO terms"""
    try:
        return(id_string.split("; "))
    except AttributeError:
        return []

def golist_flattern(nested_list):
    """Flattern a list of lists into one single big list"""
    return [item for sublist in nested_list for item in sublist]

def get_info(cluster):
    df_current = df.query('__fastGreedyCluster==@cluster')
    prot_ids = df_current['id']
    subset_edges = edges.query('source in @prot_ids | target in @prot_ids')
    identities = list(subset_edges['Identity'])

    annotation_terms = [clean_go_term(go) for go in golist_flattern(
        [extract_go_ids(id_string) for id_string in df_current['GO.IDs']]
    )]

    #get_enriched(annotation_terms)

    return [f"# Cluster_{str(cluster)}\n",
            f"Average Identity: {str(mean(identities))}\n"
            f"Enriched terms: \n",
            "## Nodes\n",
            df_current.to_markdown()+'\n',
            "## Edges\n",
            subset_edges.to_markdown()+"\n"]

if __name__ == "__main__":
    #cluster_descriptions = [get_info(cluster) for cluster in aclusters]

    for cluster in clusters:
        with open(f"clusters/cl{cluster}.md", 'w') as outfile:
            info = get_info(cluster)
            outfile.writelines(info)


#clX:
# - enriched category per cluster (versus the background)
# - output a file with clusters/cluster info
