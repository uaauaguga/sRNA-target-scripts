#!/usr/bin/env python
import argparse
from collections import defaultdict
from sklearn.metrics import roc_auc_score
import numpy as np
def main():
    parser = argparse.ArgumentParser(description='build msa of groupped sequences')
    parser.add_argument('--input', '-i', type=str, required=True, help='input scores in bed format')
    parser.add_argument('--output','-o',type=str, required=True, help="output statistics")
    args = parser.parse_args()
    leader_scores_by_genome = defaultdict(list)
    RNA_scores_by_genome = defaultdict(list)
    bg_scores_by_genome = defaultdict(list)
    with open(args.input) as f:
        for line in f:
            fields = line[:-1].split("\t")
            if len(fields) < 5:
                continue
            seq_id, score, strand = fields[0], float(fields[4]), fields[5]
            genome_id = seq_id.split(":")[0]
            seq_id = seq_id[seq_id.find(":")+1:]
            category = seq_id.split(":")[0]
            score = np.log(score/(1-score))
            if category == "leader":
                leader_scores_by_genome[genome_id].append((seq_id, score))
            elif category == "RNA":
                RNA_scores_by_genome[genome_id].append((seq_id, score))
            else:
                bg_scores_by_genome[genome_id].append(score)
    fout = open(args.output,"w")
    for genome_id in bg_scores_by_genome:
        bg_scores = bg_scores_by_genome[genome_id]
        RNA_scores = RNA_scores_by_genome[genome_id]
        leader_scores = leader_scores_by_genome[genome_id]
        #print(len(bg_scores))
        if len(bg_scores) > 10:
            mean = np.mean(bg_scores)
            std = np.std(bg_scores)
            for seq_id, score in leader_scores:
                print(genome_id, "leader", seq_id, (score - mean)/std, sep="\t", file=fout)
            for seq_id, score in RNA_scores:
                print(genome_id, "RNA", seq_id, (score - mean)/std, sep="\t",file=fout)
    fout.close()

if __name__ == "__main__":
    main()
