#!/usr/bin/env python
import numpy as np
def main():
    genome_ids = ["GCF_000016305.1","GCF_000742755.1"]
    fout = open("dataset/KP.known.sRNA.fa","w")
    with open("dataset/known-sRNAs.fa") as f:
        for line in f:  
            if line.startswith(">"):
                genome_id = line[1:].split(":")[0]
                skip = genome_id not in genome_ids                               
            if not skip:
                fout.write(line)
    fout.close()
                
            


if __name__ == "__main__":
    main() 
