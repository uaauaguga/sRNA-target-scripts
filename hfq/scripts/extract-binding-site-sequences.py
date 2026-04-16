#!/usr/bin/env python
import argparse
import pandas as pd
import subprocess
import os
from pyfaidx import Fasta
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('extract hfq associated sequences')

def main():
    parser = argparse.ArgumentParser(description='extract sequences')
    parser.add_argument('--input',  '-i', type=str, required=True, help='input interaction in bed format with pvalue at field 7')
    parser.add_argument('--genome', '-g', type=str, required=True, help='genome to extract sequence from')
    parser.add_argument('--output', '-o', type=str, help="output hfq associated sequences")
    parser.add_argument('--pvalue', '-p', type=float , default=0.001, help="pvalue cutoff")
    args = parser.parse_args()


    logger.info("Load genome ...")
    genome = Fasta(args.genome)
    genome_id = args.genome.split("/")[-1]
    genome_id = genome_id[:genome_id.rfind(".")]
    fout = open(args.output,"w") 
    logger.info("Extract sequences ...")
    with open(args.input)  as f:
        for line in f:
            fields = line[:-1].split("\t")
            pvalue = float(fields[6])
            if pvalue > args.pvalue:
                continue
            seq_id = fields[0]
            strand = fields[5]
            start, end = int(fields[1]), int(fields[2])
            p = int((start+end)/2)
            s, e = p - 50, p + 50
            s = max(0,s)
            e = min(e,len(genome[seq_id]))
            sequence = genome[seq_id][s:e]
            if strand == "-":
                sequence = sequence.reverse.complement
            sequence = str(sequence)
            print(f">{genome_id}:{seq_id}:{s}-{e}({strand})", file=fout)
            print(sequence, file=fout)             
    fout.close()
    logger.info("All done.")

if __name__ == "__main__":
    main()
