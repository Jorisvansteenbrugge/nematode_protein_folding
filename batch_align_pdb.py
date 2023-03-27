#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" This script takes two lists of pdb files and aligns them pairwise using FATCAT.
"""

__version__ = '0.1'
__author__ = "Joris J.M. van Steenbrugge"

import argparse
import logging

from itertools import product
from subprocess import Popen, STDOUT, PIPE
from joblib import Parallel, delayed


logging.basicConfig(level=logging.DEBUG)

def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--listA',
        dest='listA',
        help='list file containing the first set of pdb structures',
        required=True,
        nargs="+"
    )
    parser.add_argument(
        '--listB',
        dest='listB',
        help='list file containing the second set of pdb structures',
        required=True,
        nargs='+'
    )
    parser.add_argument(
        '--out',
        dest='out',
        help='output file',
        required=True
    )
    parser.add_argument(
        '-t',
        dest='cores',
        help='Number of parallel threads to use (default: 30). For debugging purposes use -t 1',
        default = 30,
        required=False,
        type=int
    )
    return parser.parse_args()


def parse_FATCAT_alignment(stdout):
    stdout.readline() #skip first
    scores_a = stdout.readline().decode('utf-8').split()
    scores_b = stdout.readline().decode('utf-8').split()
    
    logging.debug(scores_a)
    logging.debug(scores_b)
    out_a = [
        scores_a[1], scores_a[3], scores_a[5], scores_a[7], scores_a[9],
        scores_a[11], scores_a[13], scores_a[15], 
        ' '.join(
                [scores_a[17], scores_a[18]]
            )
        ]

    out_b = [scores_b[1], scores_b[3], scores_b[5].replace("%",""), 
             scores_b[7].replace("%","")]

    return out_a+out_b

def align_structures(prod):
    cmd = ["FATCAT", "-p1", prod[0], '-p2', prod[1], '-q'] # -q output to stdout
    logging.debug(' '.join(cmd))
    output = Popen(cmd, stdout=PIPE)


    return output.stdout

def Align(prod):
    alignment_result = align_structures(prod)
    try:
        alignment_result_parsed = parse_FATCAT_alignment(alignment_result)

        outline = list(prod)+alignment_result_parsed
    except IndexError:
        outline = [""]
    return outline

if __name__ == "__main__":
    options=get_options()

    listA = options.listA
    listB = options.listB

    logging.debug(f"List A contains: {len(listA)} values")
    logging.debug(f"List B contains: {len(listB)} values")

    FATCAT_header = [
        # First line
        "Query", "Subject", "Twists", 'ini-len', 'ini-rmsd', 'opt-equ',
        'opt-rmsd', 'chain-rmsd', 'score', 'align-len','gaps',
        # Second line
        'P-value', 'Afp-num', "Identity", "Similarity" 
    ]
    logging.debug("Using header:")
    logging.debug("\t".join(FATCAT_header))
    with open(options.out, 'w') as out:       

        alignment_lines = Parallel(
            n_jobs=options.cores)(
                delayed(Align)(prod) for prod in product(listA, listB))
        outlines = [FATCAT_header] + alignment_lines

        for line in outlines:
            out.write("\t".join(line) + '\n')    
