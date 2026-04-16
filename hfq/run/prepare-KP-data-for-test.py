#!/usr/bin/env python
import numpy as np
def main():
    genome_ids = ["GCF_000016305.1","GCF_000742755.1"]
    fout = open("dataset/test.negative.fa","w")
    with open("dataset/train.negative.fa") as f:
        for line in f:  
            if line.startswith(">"):
                genome_id = line[1:].split(":")[0]
                if np.random.rand() > 0.995:
                    skip = genome_id not in genome_ids                               
                    line = ">" + "negative:" + line[1:] 
                else:
                    skip = True
            if not skip:
                fout.write(line)
    fout.close()
                
            


if __name__ == "__main__":
    main() 
