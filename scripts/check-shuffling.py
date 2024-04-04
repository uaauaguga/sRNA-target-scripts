#!/usr/bin/env python
from ushuffle import shuffle
from scipy.stats import gumbel_r
import argparse
from multiprocessing import Pool
import subprocess
import io
from collections import defaultdict
import pickle
import numpy as np


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

def prediction(sequence_1, sequence_2, number=1, seed=7):
    cmd = ["/apps/home/lulab_jinyunfan/qhsky1/miniconda/envs/bioinfo-env/bin/IntaRNA","-q",sequence_1,"-t",sequence_2,"--outMode","C","--outNumber",str(number),"--tAcc","C","--seedBP", str(seed)]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr = subprocess.DEVNULL)
    f = io.TextIOWrapper(proc.stdout, encoding="unicode_escape")
    _ = next(f)
    records = []
    for line in f:
        if len(line.strip()) == 0:
            continue
        fields = line.strip().split(";")
        if len(fields) < 6:
            #logger.info("Some thing wrong with " + " ".join(cmd))
            continue
        t, ts, te, q, qs, qe = fields[:6]
        sequence, bp, energy = fields[-3:]
        energy = float(energy)
        ts, te, qs, qe = int(ts), int(te), int(qs), int(qe)
        records.append((qs, qe,len(sequence_1), ts, te, len(sequence_2), sequence, bp, energy))
    code =  proc.poll()
    return records

def inference5(X, params):
    X = X@params["linear_1.weight"].T + params["linear_1.bias"]
    X = np.maximum(X,0)
    X = X@params["linear_2.weight"].T + params["linear_2.bias"]
    X = np.maximum(X,0)
    X = X@params["linear_3.weight"].T + params["linear_3.bias"]
    X = np.maximum(X,0)
    X = X@params["linear_4.weight"].T + params["linear_4.bias"]
    X = np.maximum(X,0)
    return X@params["linear_5.weight"].T + params["linear_5.bias"]


pool = Pool(100)
#params5 = pickle.load(open("240423.model.pkl","rb"))
params5 = pickle.load(open("20240404.model.pkl","rb"))
def check_shuffling():    
    workers = []
    L1 = np.random.randint(50,150)
    L2 = np.random.randint(250,350)
    s1 = "".join(["ACGT"[np.random.randint(4)] for i in range(L1)])
    s2 = "".join(["ACGT"[np.random.randint(4)] for i in range(L2)])
    for i in range(10):
        sequence_1 = shuffle(s1.encode(),2).decode()
        for j in range(10):
            sequence_2 = shuffle(s2.encode(),2).decode()
            workers.append(pool.apply_async(func=prediction, args=(sequence_1, sequence_2)))
    scores = [] 
    for worker in workers:
        rs = worker.get()
        if len(rs)  == 0:
            energy = 0
        else:
            qs, qe, l1, ts, te, l2, sequence, bp, energy = rs[0]
        score = -energy/10
        scores.append(score)
    loc, scale = gumbel_r.fit(scores)
    frequency = count_frequency(sequence_1) + count_frequency(sequence_2)
    frequency = np.array(frequency).reshape(1,-1)
    pred = inference5(frequency ,params5)             
    locp, scalep = pred[0,0], pred[0,1]
    meanp = np.euler_gamma*scalep + locp
    print("#"*100)
    print("true", loc,  scale,  np.euler_gamma*scale + loc,sep="\t") 
    print("pred", locp, scalep, meanp, sep="\t")
    

def main():
    for i in range(100):
        check_shuffling()             

            


if __name__ == "__main__":
    main()
