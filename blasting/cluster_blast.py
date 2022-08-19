#!/usr/bin/env
from os import path, stat
from tempfile import NamedTemporaryFile
from jorispython import jorispython as jp
from itertools import combinations
import pandas as pd
import logging
import sys

import subprocess as sp

logging.basicConfig(level=logging.DEBUG)

clusterfile = "/home/joris/scratch/folding_comparision/cytoscape_clusters.csv"
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
                'name': 'graminicola',
                'path': "/home/joris/scratch/folding_comparision/blasting/blast_databases/meloidogyne/meloidogyne_graminicola.PRJNA411966.WBPS17.protein.fa"
            },
            {
                'name': 'hapla',
                'path': "/home/joris/scratch/folding_comparision/blasting/blast_databases/meloidogyne/meloidogyne_hapla.PRJNA29083.WBPS17.protein.fa"
            },
            {
                'name': 'incognita',
                'path': "/home/joris/scratch/folding_comparision/blasting/blast_databases/meloidogyne/meloidogyne_incognita.PRJNA340324.WBPS17.protein.fa"
            },
            {
                'name': 'chitwoodi',
                'path': "/home/joris/scratch/folding_comparision/blasting/blast_databases/meloidogyne/augustus.hints.aa"
            }
        ],
        'Gpal': [
            {
                'name': "rostochiensis",
                'path': "/home/joris/scratch/folding_comparision/blasting/blast_databases/globodera/globodera_rostochiensis.L19_PRJNA695196.WBPS17.protein.fa"
            },
            {
                'name': 'pallida',
                'path': "/home/joris/scratch/folding_comparision/blasting/blast_databases/globodera/G_pallida_D383_v.0.8.1.structural_annotation.pep.fasta"
            },
            {
                'name': "schachtii",
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
    pidentities_all = []
    for combination in combinations(clean_ids,  2):
        query_tmpfile   = NamedTemporaryFile()
        subject_tmpfile = NamedTemporaryFile()

        prep_query_fasta([combination[0]],seq_map, query_tmpfile)
        prep_query_fasta([combination[1]],seq_map, subject_tmpfile)

        output = jp.run_ncbi_blast(seq_file_name=query_tmpfile.name,
                          remote=False,
                          pairwise=True,
                          subject=subject_tmpfile.name,
                          outfmt="'6 qseqid sseqid pident qcovs evalue staxid bitscore'")

        query_tmpfile.close()
        subject_tmpfile.close()

        logging.debug(f"Blast output: {output}")

        pidentities = []
        for line in output.split("\n"):
            line = line.split()
            try:
                pidentities.append(float(line[2]))
            except:
                logging.debug(f"error in get_everage_identity line: {'-'.join(line)}")
                pass
        try:
            pidentities_all.append(sum(pidentities)/len(pidentities))
        except:
            logging.debug(f"error in division in get_everage_identity, line: {'-'.join(output)}")
            return 0
    try:
        return sum(pidentities_all)/len(pidentities_all)
    except ZeroDivisionError:
        logging.debug(f"Zero divisionerror in get_everage_identity line: {pidentities_all}")
        return None


if __name__== "__main__":

    gpal_seq_map  = jp.get_fasta_dict(gpal_sequences)
    mchit_seq_map = jp.get_fasta_dict(mchit_sequences)

    seq_maps = {'Gpal':gpal_seq_map,
                "Mchit": mchit_seq_map}


    cluster_file  = import_clusters(clusterfile)

    print("Cluster;MchitMeanIdentity;Exp_graminicola;Exp_hapla;Exp_incognita;Exp_chitwoodi;GpalMeanIdentity;Exp_rostochiensis;Exp_pallida;Exp_schachtii")

    for cluster in list(set(cluster_file["__fastGreedyCluster"])):
        logging.info(f"Starting Cl{cluster}")
        out_line = [f"Cl{cluster}"]

        for species in ["Mchit", "Gpal"]:

            tmp_file = NamedTemporaryFile()
            cluster_c = cluster_file[(cluster_file['species']== species) &
                                     (cluster_file['__fastGreedyCluster']==cluster)]

            clean_gids = [clean_cluster_gene_id(gid) for gid in list(cluster_c['name'])]

            """
            Get within species average sequence identity of the cluster
            """
            pident = get_average_identity(clean_gids, seq_maps[species])
            out_line.append(str(pident))


            """
            NCBI blast to search in other species
            """
            prep_query_fasta(
                seq_ids=clean_gids,
                seq_map=seq_maps[species],
                out_buffer=tmp_file
            )
            tmp_file.flush()
            for database in blast_databases[species]:


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
            tmp_file.close()
        print(";".join(out_line))
