#!/usr/bin/env python
import subprocess
import logging
import argparse
from scipy.stats import gumbel_r
import numpy as np
import io
import pickle
from collections import defaultdict
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('add probability')
from itertools import product
np.random.seed(666)


def count_frequency(sequence,k=2):
    counter = defaultdict(int)
    for i in range(len(sequence)-k):
        counter[sequence[i:i+k]] += 1
    frequency = np.zeros(4**k)
    idxlut = dict(zip(list("ACGT"),list(range(4))))
    for kmer in counter.keys():
        idx = 0
        skip = False
        for c in kmer:
            if c not in idxlut:
                skip = True
                break
            idx = idx*4 + idxlut[c]
        if not skip:
            frequency[idx] = counter[kmer]/10
    return list(frequency)


def inference(X, params):
    X = X@params["linear_1.weight"].T + params["linear_1.bias"]
    X = np.maximum(X,0)
    X = X@params["linear_2.weight"].T + params["linear_2.bias"]
    X = np.maximum(X,0)
    X = X@params["linear_3.weight"].T + params["linear_3.bias"]
    X = np.maximum(X,0)
    X = X@params["linear_4.weight"].T + params["linear_4.bias"]
    X = np.maximum(X,0)
    return X@params["linear_5.weight"].T + params["linear_5.bias"]


def main():
    parser = argparse.ArgumentParser(description='calculate interaction probability')
    parser.add_argument('--input',  '-i', type=str, required=True, help='Input interaction energy')
    parser.add_argument('--output','-o', type=str , help="where to save output")
    parser.add_argument('--model', '-m', type=str, default = "20240404.model.pkl", help='model to use')
    parser.add_argument('--word-size', '-k', type=int, default = 2, help='word size to use')
    args = parser.parse_args()
    params = pickle.load(open(args.model,"rb")) 
    frequencies = []
    energies = []
    cached_lines = []
    fout = open(args.output,"w")
    n = 0
    with open(args.input) as f:
        for line in f:
            sequence_1, sequence_2, energy = line.strip().split("\t")[:3]
            cached_lines.append("\t".join([sequence_1, sequence_2, energy]))
            energies.append(float(energy))
            frequency = count_frequency(sequence_1,args.word_size) + count_frequency(sequence_2,args.word_size) 
            frequencies.append(np.array(frequency))
            n += 1
            if len(frequencies) == 4096:
                assert len(cached_lines) == 4096
                X = np.array(frequencies)
                X = inference(X, params) 
                for i, line in enumerate(cached_lines):
                    loc, scale = X[i,0], X[i,1]
                    energy = energies[i]
                    prob = gumbel_r.cdf(-energy/10, loc=loc,scale=scale)
                    print(line, prob, sep="\t", file=fout)
                energies = []
                frequencies = []
                cached_lines = []
            if n%10000 == 0:
                logger.info(f"{round(n/1000)} K pairs processed .")
    if len(frequencies) > 0:
        X = np.array(frequencies)
        X = inference(X, params)
        for i, line in enumerate(cached_lines):
            loc, scale = X[i,0], X[i,1]
            energy = energies[i]
            prob = gumbel_r.cdf(-energy/10, loc=loc, scale=scale)
            print(line, prob, sep="\t", file=fout)
    fout.close()
    logger.info("all done .")

if __name__ == "__main__":
    main()
