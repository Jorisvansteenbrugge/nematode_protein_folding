#!/usr/bin/env python3
import pandas as pd
import subprocess as sp
from tempfile import NamedTemporaryFile


def get_genes_cluster(cluster_name, network_table, species_marker):
    sub = network_table[network_table["__fastGreedyCluster"] == cluster_name]
    return list(sub[sub["species"] == species_marker]["shared name"])

def clean_gene_names(cluster, network_table, species_marker):
    genes = get_genes_cluster(cluster, network_table, species_marker)
    return [gene.replace(".pdb", ".t1") for gene in genes]

def _has_interGO(value):
    return value not in ("no IPS match", "no GO terms")

def _has_description(value):
    if "-NA-" in value or 'unnamed protein' in value or 'hypothetical protein' in value:
        return False
    else:
        return True

def is_annotated(row, t):
    description = _has_description(row["Description"])
    interGO = _has_interGO(row["InterPro GO IDs"])
    return False in (description, interGO)

def check_annotations(fun_anno, spprefix, annotation_table):
    out_table = {}
    annotated = list(
        fun_anno.apply(is_annotated,
                            axis=1,
                            t=annotation_table))
    seqids = list(fun_anno["SeqName"])


    unannotated = [not x for x in annotated]
    unannotated_gene_ids = list(fun_anno[unannotated]["SeqName"])


    if False in annotated:
        if True in annotated: #Knowable
            return [(seqids[i], annotated[i]) for i in range(len(annotated)) ]
        else: return False #unkown
    else: return False # all known not interesting

def get_matches(gid: str):
    return [x.replace(".pdb",".t1") for x in list(edge_table.query("target==@gid")["source"])]

def export_clusters(cluster: str, species: str, checked: list):
    with open(f"Cl{cluster}_knowable_{species}.csv", 'w') as outfile:
        for line in checked:
            outfile.write(";".join([str(x) for x in line])+'\n')
def export_knowable_table(cluster: str, species: str, checked: list, handler):
    for line in checked:
        if line[1]:
            gid = line[0]
            gid_pdb = gid.replace(".t1",".pdb")
            handler.write(f"{cluster}\t{species}\t{gid}\t{gid}\n")
            for match in get_matches(gid_pdb):
                handler.write(f"{cluster}\t{species}\t{gid}\t{match}\n")

def extract_fastas():
    id_cmd = 'cat knowable_cluster_table.csv | cut -f 4 | tail -n +2'
    ids = sp.check_output(id_cmd, shell = True).decode().split('\n')

    cmd = [
        'get_fasta_regions.py',
        '--fasta',
        '/mnt/nemahomes/steen176/Genomes/Annotations/G_pallida_D383_v.0.8.1.structural_annotation.pep.fasta',
        "/mnt/dataHDD/scratch/chitwoodi_genomes/annotation/braker/augustus.hints.aa",
        '--ids' ] + ids + [">", "knowable_cluster_table_sequences.fasta"]

    sp.call(" ".join(cmd), shell = True)

def hmmscan():
    cmd = ["bash","hmmscan_script.sh", "knowable_cluster_table_sequences.fasta"]
    sp.call(cmd)


def check_evalue(line):
    return True if float(line[6]) < 0.1 else False

def import_hmmscan():
    domain_map = {}
    with open("knowable_cluster_table_sequences.domtblout") as domtbl:
        for line in domtbl:
            if line.startswith("#"):
                continue

            line = line.strip().split()

            if not check_evalue(line):
                continue

            query = line[3]
            desc= " ".join(line[22:])
            target_n_desc = (line[0],desc)

            if query not in domain_map:
                domain_map[query] = [target_n_desc]
            else:
                domain_map[query].append(target_n_desc)
    return domain_map


def create_knowable_annotation_cluster_table(domain_map,outf):
    with open(outf,'w') as outfile:
        with open("knowable_cluster_table.csv") as table:
            header = table.readline().strip().split("\t") + ["Annotation", "Pfam Domains"]
            outfile.write("\t".join(header) + "\n")

            for line in table:
                line = line.strip().split()
                target = line[3]

                if target.startswith('Gpal_D383'):
                    gene_annotation = list(gp_gene_annotations.query("SeqName==@target")["Description"])[0]
                else:
                    gene_annotation = list(mc_gene_annotations.query("SeqName==@target")["Description"])[0]


                try:
                    domains = ';'.join([','.join(x) for x in domain_map[target]])
                except KeyError:
                    domains = 'NA'
                outfile.write("\t".join(line)+f"\t{gene_annotation}\t{domains}\n")


    
network_table = pd.read_csv("../cytoscape_clusters_new.csv",
                            sep=',')

gp_gene_annotations = pd.read_csv("/mnt/nemahomes/steen176/Genomes/Annotations/G_pallida_annotations.csv",
                                   sep="\t")

mc_gene_annotations = pd.read_csv("/mnt/dataHDD/scratch/chitwoodi_genomes/functional_annotation2.csv",
                                  sep='\t')
edge_table = pd.read_csv("../cytoscape_edges.csv",
                         sep = ';')

with open("knowable_cluster_table.csv",'w') as knowable_table:
    knowable_table.write("Cluster\tSpecies\tUnannotated_gene\tTarget\n")
    clusters = list(set(network_table["__fastGreedyCluster"]))
    for cluster in clusters:
        gp_genes = clean_gene_names(cluster, network_table, "Gpal")
        gp_functional = gp_gene_annotations.query("SeqName.isin(@gp_genes)")

        mc_genes = clean_gene_names(cluster, network_table, "Mchit")
        mc_functional = mc_gene_annotations.query("SeqName.isin(@mc_genes)")

        gp_checked = check_annotations(gp_functional, "Gpal", gp_gene_annotations)
        mc_checked = check_annotations(mc_functional, "Mchit", mc_gene_annotations)

        if isinstance(gp_checked, list):
            export_clusters(str(cluster),'Gpal', gp_checked)
            export_knowable_table(str(cluster),"Gpal", gp_checked,knowable_table)
        if mc_checked != False:
            export_clusters(str(cluster), "Mchit", mc_checked)

extract_fastas()
hmmscan()
domain_map = import_hmmscan()
create_knowable_annotation_cluster_table(domain_map, outf="knowable_cluster_table_annnotated.csv")