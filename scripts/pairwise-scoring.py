#!/usr/bin/env python
import subprocess
import logging
import argparse
from scipy.stats import gumbel_r
import numpy as np
from pyfaidx import Fasta
import io
import pickle
from collections import defaultdict
from multiprocessing import Pool
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('target prediction')
from itertools import product
np.random.seed(666)

"""
def count_frequency(sequence,k=2):
    counter = defaultdict(int)
    for i in range(len(sequence)-k):
        counter[sequence[i:i+k]] += 1
    frequency = []
    for kmer in list(product(*k*["ACGT"])):
        kmer = "".join(kmer)
        frequency.append(counter[kmer]/10)
    return frequency
"""

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


def prediction(sequence_1, sequence_2, number, seed):
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

def load_fasta(path):
    sequences = defaultdict(dict)
    with open(path) as f:
       for line in f:
           if line.startswith(">"):
               seq_id = line[1:].strip().split(" ")[0]
               gene_id, genome_id = seq_id.split("--")
               ref_genome_id = gene_id.split(":")[0]
               target_genome_id = genome_id.split(":")[0]
               sequences[(ref_genome_id, gene_id)][target_genome_id] = ""
           else:
               sequences[(ref_genome_id, gene_id)][target_genome_id] += line.strip()
    return sequences


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
    parser = argparse.ArgumentParser(description='predict RNA RNA interaction with intaRNA')
    parser.add_argument('--srnas',  '-rs', type=str, required=True, help='sRNA sequences')
    parser.add_argument('--targets',  '-ts', type=str, required=True, help='Target sequences')
    parser.add_argument('--output','-o', type=str , help="where to save output")
    parser.add_argument('--number','-n', type=int , default=1, help="number of prediction per sequence pair")
    parser.add_argument('--jobs','-j', type=int , default=128, help="number of process to run")
    parser.add_argument('--seed','-s', type=int , default=7, help="seed length to use") 
    parser.add_argument('--model', '-m', type=str, default = "240423.model.pkl", help='model to use')
    parser.add_argument('--word-size', '-k', type=int, default = 2, help='word size to use')
    args = parser.parse_args()
    
    fout = open(args.output,"w")
    logger.info("load sequences ...")
    sRNAs = load_fasta(args.srnas)
    targets = load_fasta(args.targets) 
    logger.info(f"run intaRNA with {args.jobs} workers ...")
    pool = Pool(args.jobs)
    n_too_long = 0
    n_no_prediction = 0
    n_total = 0
    n_too_short = 0
    n_few_reads = 0
    logger.info("load weights for background distribution modeling ...")
    params = pickle.load(open(args.model,"rb"))


    print("sRNA id","sstart","send","target id","tstart","tend","energy","pvalue","# intersection","# union","loc","scale",file=fout,sep="\t")    
    for qref_genome_id, sRNA_id in sRNAs:
        for tref_genome_id, target_id in list(targets.keys()):
            if qref_genome_id != tref_genome_id:
                continue
            ref_genome_id = qref_genome_id 
            workers = []
            frequencies = []
            sRNA_genome_ids = list(sRNAs[(ref_genome_id,sRNA_id)].keys())
            target_genome_ids = list(targets[(ref_genome_id,target_id)].keys())
            paired_genome_ids = np.intersect1d(sRNA_genome_ids,target_genome_ids)
            if len(paired_genome_ids) <= 6:
                continue
            n_sRNA_unique_genome = len(sRNA_genome_ids) - len(paired_genome_ids)
            n_target_unique_genome = len(target_genome_ids) - len(paired_genome_ids)
            n_union_genome = len(paired_genome_ids) + n_sRNA_unique_genome + n_target_unique_genome
            if len(paired_genome_ids) > 4096:
                np.random.shuffle(paired_genome_ids)
                paired_genome_ids = paired_genome_ids[:4096]
            paired_genome_ids = sorted(list(paired_genome_ids))
            for genome_id in paired_genome_ids:
                sequence_1 = sRNAs[(ref_genome_id,sRNA_id)][genome_id]
                sequence_2 = targets[(ref_genome_id,target_id)][genome_id]
                workers.append((pool.apply_async(func=prediction, args=(sequence_1, sequence_2, args.number,args.seed)),sRNA_id, target_id, genome_id, len(paired_genome_ids), n_union_genome))
                frequency = count_frequency(sequence_1,args.word_size) + count_frequency(sequence_2,args.word_size) 
                frequencies.append(np.array(frequency))
            if len(frequencies) == 0:
                continue
            X = np.array(frequencies)
            X = inference(X, params) 
            i = 0
            logger.info(f"{sRNA_id} {target_id}: {len(workers)} interactions to process .")
            n_no_prediction = 0
            for worker, sRNA_id, target_id, genome_id, n_paired, n_union in workers:
                seq_id_1 = sRNA_id + ":" + genome_id
                seq_id_2 = target_id + ":" + genome_id
                rs = worker.get()
                loc, scale = X[i,0], X[i,1]
                i += 1
                for r in rs:
                    qs, qe, l1, ts, te, l2, sequence, bp, energy = r
                    pvalue = 1 - gumbel_r.cdf(-energy/10, loc=loc,scale=scale) 
                    print(seq_id_1, qs-1, qe, seq_id_2, ts-1, te, energy, pvalue, n_paired, n_union,  loc, scale, sep="\t", file=fout)
                if len(rs) == 0:
                    n_no_prediction += 1
                    print(seq_id_1,"-1","-1",seq_id_2,"-1","-1",0,1, n_paired, n_union, sep="\t",file=fout)            
            fout.flush()
    logger.info(f"{n_no_prediction} have no prediction.")
    fout.close()
    logger.info("all done .")

if __name__ == "__main__":
    main()
