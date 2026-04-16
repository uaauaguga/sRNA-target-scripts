#!/usr/bin/env python
from collections import defaultdict
import numpy as np
def main():
    scores = defaultdict(list)
    with open("Rfam.seed.bed") as f:
        for line in f:
            fields = line.strip().split("\t")
            score = float(fields[4])
            rfam_id, rfam_name = fields[0].split(":")[:2]
            rfam_acc = rfam_id + "-" + rfam_name
            scores[rfam_acc].append(score)
    rfam_ids = open("rfam-bacteria.txt").read().strip().split("\n")
    for rfam_acc in rfam_ids:#scores:
        print(rfam_acc, np.median(scores[rfam_acc]),sep="\t") 

if __name__ == "__main__":
    main()
