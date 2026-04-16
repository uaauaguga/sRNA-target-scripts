#!/usr/bin/env python
import os
import torch
from dataset import onehot
from torch.functional import F
import numpy as np
from model import CNNClassifier
import argparse
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('scanning sequence')
from tqdm import tqdm
from scipy.signal import convolve

rc_lut = {"A":"T","C":"G","G":"C","T":"A"}


def get_rc(sequence):
    rcs = []
    for c in sequence:
        c = rc_lut.get(c,'N')
        rcs.append(c)
    return "".join(rcs[::-1])

def main():
    parser = argparse.ArgumentParser(description='subsequence scoring')
    parser.add_argument('--fasta','-f',type=str, required=True, help="Input sequence")
    parser.add_argument('--batch-size','-bs',type=int,default=1024,help="Batch size for scanning")
    parser.add_argument('--device','-d',default="cuda:0",choices=["cuda:0","cuda:1","cuda:2","cuda:3","cpu"],help="Device to run the model")
    parser.add_argument('--model','-m', default = "models.b10.c64.with.targets.da0.5/36.pt", type=str,help="Where to load the model parameters")
    parser.add_argument('--output','-o',required=True,type=str,help="Where to save the predicted probabilities")
    parser.add_argument('--n-channels', '-nc', type=int, default=64, help="number of channels to use")
    parser.add_argument('--n-blocks', '-nb', type=int, default=10, help="number of blocks in the model")
    parser.add_argument('--kernel-size', '-k', type=int, default=5, help="convolution kernel size in the res-block")
    parser.add_argument('--length', '-l', type=int, default=100, help="length of an instance")
    parser.add_argument('--reverse-complementary', '-rc', action = "store_true", help="scoring both strand")
    parser.add_argument('--stride', '-s', type=int, default = 1, help="stride size for scanning")
    parser.add_argument('--cutoff', '-c', type=float, default = 0, help="cutoff")
    parser.add_argument('--offset', type=int, default = 0, help="offset relative to predicted ends")
    args = parser.parse_args()

    from pyfaidx import Fasta
    logger.info("load sequence ...")
    fasta = Fasta(args.fasta)
    seq_ids, sequences = [], []
    for seq_id in fasta.keys():
        #logger.info(f"load {seq_id} ...")
        seq_ids.append(seq_id)
        sequences.append(str(fasta[seq_id][::]))
    logger.info(f"load model weights from {args.model} ...")
    state_dict = torch.load(args.model, map_location = args.device)
    if 'daOut.weight' in state_dict:
        n_domains = state_dict['daOut.weight'].shape[0]
    else:
        n_domains = None
    model = CNNClassifier(n_domains = n_domains, n_channels=args.n_channels, kernel_size = args.kernel_size, n_blocks=args.n_blocks)
    model = model.to(args.device)
    model = model.eval()
    model.load_state_dict(state_dict) 
    model.daOut = None

    batched_instances = []
    batched_positions = []

    positions_by_sequence = {}
    scores_by_sequence = {}

    logger.info("run prediction ...")
    n = 0
    for seq_id, sequence in zip(seq_ids, sequences):
        if (n+1)%10000 == 0:
            logger.info(f"{round(n/1000)} K sequence processed")
        L = len(sequence)
        n += 1
        if L < args.length:
            sequence = sequence + "N"*(args.length-len(sequence))
            L = args.length
        sequence = sequence.upper().replace("U","T")
        if args.reverse_complementary:
            rc_sequence = get_rc(sequence)
            rc_sequence = onehot(rc_sequence)
        sequence = onehot(sequence)
        p = 0
        while p+args.length <= L:
            instance = sequence[:,:,p:p+args.length]
            batched_instances.append(instance)
            batched_positions.append((seq_id, p,"+"))
            if args.reverse_complementary:
                instance = rc_sequence[:,:,L-p-args.length:L-p]
                batched_instances.append(instance)
                batched_positions.append((seq_id, p,"-"))
            p += args.stride
            if len(batched_instances) == args.batch_size:
                X = torch.cat(batched_instances).to(args.device)
                y = model(X)                
                scores = torch.softmax(y,1)[:,1].detach().cpu().numpy()
                for i, (seq_id, position, strand) in enumerate(batched_positions):
                    score = scores[i]
                    #if score < args.cutoff:
                    #    continue
                    if (seq_id, strand) not in scores_by_sequence:
                        scores_by_sequence[(seq_id, strand)] = []
                        positions_by_sequence[(seq_id, strand)] = []
                    scores_by_sequence[(seq_id, strand)].append(score)
                    positions_by_sequence[(seq_id, strand)].append(position)
                batched_instances = []
                batched_positions = []
    if len(batched_instances) > 0:
        X = torch.cat(batched_instances).to(args.device)
        y = model(X)
        scores = torch.softmax(y,1)[:,1].detach().cpu().numpy()        
        for i, (seq_id, position, strand) in enumerate(batched_positions):
            score = scores[i]
            #if score < args.cutoff:
            #    continue
            if (seq_id, strand) not in scores_by_sequence:
                scores_by_sequence[(seq_id, strand)] = []
                positions_by_sequence[(seq_id, strand)] = []
            scores_by_sequence[(seq_id, strand)].append(score)
            positions_by_sequence[(seq_id, strand)].append(position)

    
    logger.info("processing the predictions ...")
    logger.info(f"results will be saved to {args.output} ...")
    fout = open(args.output,"w")

    rfam_ids = set()
    for seq_id, strand in scores_by_sequence:
        rfam_ids.add(seq_id.split(":")[0])
        scores = scores_by_sequence[(seq_id, strand)]
        positions = positions_by_sequence[(seq_id, strand)]
        for position,score in zip(positions,scores):
            if strand == "+":
                print(seq_id, position + args.offset, position + args.length + args.offset, ".", score, strand, sep="\t", file=fout)
            else:
                print(seq_id, position + args.offset, position + args.offset + args.length, ".", score, strand, sep="\t", file=fout)
    print(len(rfam_ids))
    fout.close()
    logger.info("all done .")

                
if __name__ == "__main__":
    main()
