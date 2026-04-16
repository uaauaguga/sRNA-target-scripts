#!/usr/bin/env python
import argparse
from collections import defaultdict
from sklearn.metrics import roc_auc_score, roc_curve
import numpy as np
import pandas as pd
def main():
    parser = argparse.ArgumentParser(description='build msa of groupped sequences')
    parser.add_argument('--input', '-i', type=str, required=True, help='input scores in bed format')
    parser.add_argument('--output','-o',type=str, required=True, help="output statistics")
    args = parser.parse_args()
    scores_by_categories = defaultdict(list)
    labels = []
    scores = []
    with open(args.input) as f:
        for line in f:
            fields = line[:-1].split("\t")
            seq_id, score, strand = fields[0], float(fields[4]), fields[5]
            label = "negative" not in seq_id
            labels.append(int(label))
            scores.append(score)
        AUROC = roc_auc_score(labels,scores)
    print(AUROC)
    tpr, fpr, threshold = roc_curve(labels,scores)
    ROC_curve = {"tpr":tpr,"fpr":fpr,"threshold":threshold}
    pd.DataFrame(ROC_curve).to_csv(args.output, index=False,sep=",")

if __name__ == "__main__":
    main()
