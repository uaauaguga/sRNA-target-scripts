#!/usr/bin/env python
import numpy as np
def main():
    genome_ids = ["GCF_000016305.1","GCF_000742755.1"]
    ftrain = open("dataset/RIL-seq.targets.wo.KP.fa","w")
    ftest = open("dataset/RIL-seq.targets.KP.fa","w")
    with open("../hfq-binding-prediction/dataset/target-sequences.99.9.fa") as f:
        for line in f:  
            if line.startswith(">"):
                genome_id = line[1:].split(":")[0]
                if genome_id not in genome_ids:
                    fout = ftrain
                else:
                    fout = ftest
            fout.write(line)
    ftrain.close()
    ftest.close()
                
            


if __name__ == "__main__":
    main() 
