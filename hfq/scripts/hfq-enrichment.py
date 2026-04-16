#!/usr/bin/env python
import argparse
from collections import defaultdict
from sklearn.metrics import roc_auc_score
import numpy as np
def main():
    parser = argparse.ArgumentParser(description='build msa of groupped sequences')
    parser.add_argument('--input', '-i', type=str, required=True, help='input scores in bed format')
    parser.add_argument('--output','-o',type=str, required=True, help="output statistics")
    parser.add_argument('--n-bins','-nb',type=int, default=1000, help="number of bins")
    parser.add_argument('--n-levels','-nl',type=int, default=20, help="number of levels")
    args = parser.parse_args()
    scores_by_categories = defaultdict(list)
    with open(args.input) as f:
        for line in f:
            fields = line[:-1].split("\t")
            seq_id, score, strand = fields[0], float(fields[4]), fields[5]
            seq_id = seq_id[seq_id.find(":")+1:]
            category = seq_id.split(":")[0]
            scores_by_categories[category].append(score)
            if category == "RNA":
                length = int(seq_id.split(":")[-1])
                if length < 40:
                    scores_by_categories["RNA.lt50"].append(score)
                else:
                    scores_by_categories["RNA.ge50"].append(score)
    for k in ["RNA","RNA.lt50","RNA.ge50"]:
        y_true = len(scores_by_categories["random"])*[0] + len(scores_by_categories[k])*[1]
        y_pred = scores_by_categories["random"] + scores_by_categories[k]
        AUROC = roc_auc_score(y_true,y_pred)
        print(f"AUROC score of {k} is: {AUROC}.")
    score2level = {}
    n_levels = args.n_levels
    n_bins = args.n_bins
    for i,score in enumerate(sorted(scores_by_categories['random'])[::-1]):
        score = round(n_bins*score)
        level = int(i*n_levels/len(scores_by_categories['random']))
        score2level[score] = level

    score2level[0] = n_levels
    for i in range(1,n_bins):
        if i not in score2level:
            score2level[i] = score2level[i-1]
    leader_score_distribution = np.zeros(n_levels+1)
    for score in scores_by_categories['leader']:
        score = int(score*n_bins)
        level = score2level[score]
        leader_score_distribution[level] += 1
    leader_score_distribution = leader_score_distribution/leader_score_distribution.sum()
    RNA_score_distribution = np.zeros(n_levels+1)
    for score in scores_by_categories['RNA']:
        score = int(score*n_bins)
        level = score2level[score]
        RNA_score_distribution[level] += 1
    RNA_score_distribution = RNA_score_distribution/RNA_score_distribution.sum()
    random_score_distribution = np.zeros(n_levels+1)
    for score in scores_by_categories['random']:
        score = int(score*n_bins)
        level = score2level[score]
        random_score_distribution[level] += 1
    random_score_distribution = random_score_distribution/random_score_distribution.sum()

    RNA_enrichment = (RNA_score_distribution+0.01)/(random_score_distribution+0.01)
    leader_enrichment = (leader_score_distribution+0.01)/(random_score_distribution+0.01)

    fout = open(args.output,"w")
    print("category","RNA","leader",sep="\t",file=fout)
    for i in range(n_levels):
        print(i, RNA_enrichment[i], leader_enrichment[i], sep="\t",file=fout)
        print(i, RNA_enrichment[i], leader_enrichment[i], sep="\t")

    #for score in sorted(list(score2level.keys())):
    #    print(score, score2level[score], sep="\t")

if __name__ == "__main__":
    main()
