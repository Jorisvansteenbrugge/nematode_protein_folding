#!/usr/bin/env python3

from tempfile import NamedTemporaryFile
from jorispython import jorispython as jp
import numpy as np
import logging

logging.basicConfig(level=logging.INFO)


gpal_seqs = jp.get_fasta_dict(
    "/mnt/nemahomes/steen176/Genomes/Annotations/G_pallida_D383_v.0.8.1.structural_annotation.pep_rbpfix.fasta"
)
mchit_seqs = jp.get_fasta_dict(
    "/mnt/dataHDD/scratch/chitwoodi_genomes/annotation/braker/augustus.hints.aa"
)

def to_tmp_file(sid,map, outf):
    sid = sid.replace(".pdb",".t1")
    seq = map[sid]
    fasta_entry = f">{sid}\n{seq}\n"
    outf.write(fasta_entry.encode())
    outf.flush()

with open("folding_fatcat_output_rbp_vs_all.tsv") as infile:
    header = infile.readline().strip().split('\t')
    header.append("global_identity")
    print("\t".join(header))

    c = 1

    for line in infile:
        if c%100 == 0:
            logging.info(f"Done: {c}")

        line = line.strip().split("\t")

        query = line[0]
        subject = line[1]


        query_tmp = NamedTemporaryFile()
        subject_tmp = NamedTemporaryFile()



        to_tmp_file(query, mchit_seqs, query_tmp)
        to_tmp_file(subject, gpal_seqs, subject_tmp)

        logging.debug(query_tmp.name)


        needle_out = jp.run_pairwise_needle(query_tmp.name, subject_tmp.name)
        pidentity = jp.get_identity_from_needle_output(needle_out)


        line.append(pidentity)

        print("\t".join(line))
        c+=1
