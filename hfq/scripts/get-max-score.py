#!/usr/bin/env python
from collections import defaultdict
import numpy as np
def main():
    scores = defaultdict(list)
    #with open("test.scores.bed") as f:
    #with open("test.scores.5.bed") as f:
    #with open("trans.acting.sRNA.5.bed") as f:
    #with open("test.5.bed") as f:
    with open("2.5.bed") as f:
        for line in f:
            fields = line.strip().split("\t")
            start = int(fields[1])
            if start > 10:
                continue
            seq_id = fields[0]
            score = float(fields[4])
            scores[seq_id].append(score)
    for seq_id in scores:
        print(seq_id, np.max(scores[seq_id]),sep="\t") 

if __name__ == "__main__":
    main()
