#!/usr/bin/env python
import numpy as np
np.random.seed(666)
def main():
    genome_ids = open("proteobacteria-sampled-genomes.txt").read().strip().split("\n")
    genome_ids = set(genome_ids)
    fout = open("dataset/inference/GTDB.proteobacteria.fa","w")
    for i in range(623):
        print(i)
        chunk_id = str(i).zfill(4)       
        path = f"dataset/inference/GTDB/{chunk_id}.fa"
        with open(path) as f:
            for header in f:
                sequence = next(f)
                genome_id = header[1:].strip().split(":")[0]
                if genome_id in genome_ids:
                    if np.random.rand() < 0.1:
                        fout.write(header)
                        fout.write(sequence)
    fout.close()
if __name__ == "__main__":
    main()
