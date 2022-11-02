#!/usr/bin/env
from os import path, stat
from tempfile import NamedTemporaryFile
from jorispython import jorispython as jp
from itertools import combinations
import pandas as pd
import logging
import sys
import numpy as np
import subprocess as sp

logging.basicConfig(level=logging.INFO)

clusterfile = "/home/joris/scratch/folding_comparision/cytoscape_clusters_new.csv"
"""
Index(['__fastGreedyCluster', 'Description', 'Enzyme.Names', 'GO.IDs', 'id',
       'InterPro.GO.Names', 'InterPro.IDs', 'name', 'selected', 'shared name',
       'species'],
      dtype='object')
"""

gpal_sequences = "/mnt/nemahomes/steen176/Genomes/Annotations/G_pallida_D383_v.0.8.1.structural_annotation.pep.fasta"  #
"""
>Gpal_D383_g00001.t1
MAVVFILMLLSVATLEAQLRNEKDEFFDLMTAIGEYFGGEFNDAEKALELDLNGANSFLE
"""

mchit_sequences = "/mnt/dataHDD/scratch/chitwoodi_genomes/annotation/braker/augustus.hints.aa"
"""
>g6761.t1
MSRISTAVSSHFDIQSNAIFEPKTALHRQFEKIKSSKETFKNNMQIEPHKVKELELILTE
"""

blast_databases = {
        'Mchit': [
            {
                'name': 'M. graminicola',
                'path': "/home/joris/scratch/folding_comparision/blasting/blast_databases/meloidogyne/meloidogyne_graminicola.PRJNA411966.WBPS17.protein.fa"
            },
            {
                'name': 'M. hapla',
                'path': "/home/joris/scratch/folding_comparision/blasting/blast_databases/meloidogyne/meloidogyne_hapla.PRJNA29083.WBPS17.protein.fa"
            },
            {
                'name': 'M. incognita',
                'path': "/home/joris/scratch/folding_comparision/blasting/blast_databases/meloidogyne/meloidogyne_incognita.PRJNA340324.WBPS17.protein.fa"
            },
            {
                'name': 'M. chitwoodi',
                'path': "/home/joris/scratch/folding_comparision/blasting/blast_databases/meloidogyne/augustus.hints.aa"
            }
        ],
        'Gpal': [
            {
                'name': "G. rostochiensis",
                'path': "/home/joris/scratch/folding_comparision/blasting/blast_databases/globodera/globodera_rostochiensis.L19_PRJNA695196.WBPS17.protein.fa"
            },
            {
                'name': 'G. pallida',
                'path': "/home/joris/scratch/folding_comparision/blasting/blast_databases/globodera/G_pallida_D383_v.0.8.1.structural_annotation.pep.fasta"
            },
            {
                'name': "H. schachtii",
                'path': "/home/joris/scratch/folding_comparision/blasting/blast_databases/globodera/H_schachtii_IRS_PROT_named.fasta"
            }
        ]
    }

def import_clusters(cluster_file_name, sep=','):
    return pd.read_csv(cluster_file_name, sep=sep)

def prep_query_fasta(seq_ids, seq_map, out_buffer):
    logging.debug(f'prep_query_fasta input ids was: {",".join(seq_ids)}')
    out = [f">{seq_id}\n{seq_map[seq_id]}" for seq_id in seq_ids]
    logging.debug(f'prep_query_fasta output list looks like this: {";".join(out)}')
    out_buffer.write("\n".join(out).encode())
    out_buffer.flush()

def clean_cluster_gene_id(cluster_gene_id):
    return cluster_gene_id.replace('.pdb','.t1')

def has_too_many_sequences(gene_ids, threshold=25):
    return True if len(gene_ids) > 25 else False

def filter_blast_output(dataframe, evalue_thresh=0.05, pident_thresh=70,qcov_thresh=70):
    """Criteria will be set to:
        evalue < 0.05
        pidentity >= 35%
        qcov_thresh >= 70%
    """
    dataframe['pident'] = pd.to_numeric(dataframe['pident'])
    dataframe['qcovs']  = pd.to_numeric(dataframe['qcovs'])
    dataframe['evalue'] = pd.to_numeric(dataframe['evalue'])

    filtered =dataframe[
        (dataframe['pident'] >= pident_thresh) &
        (dataframe['qcovs']  >= qcov_thresh) &
        (dataframe['evalue'] <= evalue_thresh)
    ]
    return filtered

def get_average_identity(clean_ids, seq_map):
    if len(clean_ids) == 1:
        return np.nan

    pidentities_all = []
    for combination in combinations(clean_ids,  2):
        query_tmpfile   = NamedTemporaryFile()
        subject_tmpfile = NamedTemporaryFile()

        prep_query_fasta([combination[0]],seq_map, query_tmpfile)
        prep_query_fasta([combination[1]],seq_map, subject_tmpfile)

        raw_needle_output = jp.run_pairwise_needle(query_tmpfile.name,subject_tmpfile.name)
        pairwise_identity = jp.get_identity_from_needle_output(raw_needle_output)

        query_tmpfile.close()
        subject_tmpfile.close()

        logging.info(f"Pairwise Identity: {pairwise_identity}")

        pidentities_all.append(float(pairwise_identity))


    try:
        return (sum(pidentities_all)/len(pidentities_all))
    except ZeroDivisionError:
        #print(pidentities_all)
        loggin.ERROR("ERROR")
        sys.exit()



if __name__== "__main__":

    gpal_seq_map  = jp.get_fasta_dict(gpal_sequences)
    mchit_seq_map = jp.get_fasta_dict(mchit_sequences)

    seq_maps = {'Gpal':gpal_seq_map,
                "Mchit": mchit_seq_map}


    cluster_file  = import_clusters(clusterfile)

#    print("Cluster;MchitMeanIdentity;Exp_graminicola;Exp_hapla;Exp_incognita;Exp_chitwoodi;GpalMeanIdentity;Exp_rostochiensis;Exp_pallida;Exp_schachtii")
    print("Cluster;Species;Mean_identity;Exp_species;Exp_value")

    for cluster in list(set(cluster_file["__fastGreedyCluster"])):
        logging.info(f"Starting Cl{cluster}")


        for species in ["Mchit", "Gpal"]:


            tmp_file = NamedTemporaryFile()

            cluster_c = cluster_file[(cluster_file['species']== species) &
                                     (cluster_file['__fastGreedyCluster']==cluster)]

            clean_gids = [clean_cluster_gene_id(gid) for gid in list(cluster_c['name'])]
            logging.debug(clean_gids)
            """
            Get within species average sequence identity of the cluster
            """
            pident = get_average_identity(clean_gids, seq_maps[species])
            logging.info(f"average pairwise identity Cl{cluster} for {species}: {pident}")
            #out_line.append(str(pident))



            """
            NCBI blast to search in other species
            """
            logging.debug(f"Found gene ids: {','.join(clean_gids)}")
            prep_query_fasta(
                seq_ids=clean_gids,
                seq_map=seq_maps[species],
                out_buffer=tmp_file
            )
            tmp_file.flush()
            for database in blast_databases[species]:
                out_line = [f"Cl{cluster}", species, str(pident), database['name']]

                raw_blastout = jp.run_ncbi_blast(
                    seq_file_name=tmp_file.name,

                    remote=False,
                    db=database['path']
                    )
                logging.debug(f"Blast output Cl{cluster} against {database['name']}: {raw_blastout}")
                parsed_blastout = [line.split('\t') for line in raw_blastout.split('\n')]

                try:
                    blast_df = pd.DataFrame(parsed_blastout,
                                        columns = ['qseqid','sseqid','pident','qcovs','evalue','staxid']
                                        )
                    filtered_blast_output = filter_blast_output(blast_df)
                    unique_target_seqids = list(set(filtered_blast_output['sseqid']))
                    logging.info(f"Cl{cluster} matched with the following {database['name']} gene ids {unique_target_seqids}")

                    out_line.append(str(len(unique_target_seqids)))
                except ValueError:
                    out_line.append('0')

                """
                Output statistics to csv
                """
                print(";".join(out_line))

            tmp_file.close()
