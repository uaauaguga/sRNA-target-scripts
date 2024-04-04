#!/usr/bin/env python
from itertools import product
import torch
from torch.utils.data import Dataset, DataLoader, SequentialSampler,RandomSampler
from torch.functional import F
from torch import nn
from collections import defaultdict
from torch.optim import Adam
import numpy as np
import sys
import argparse
import os
import re
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] %(message)s')
logger = logging.getLogger("model training")


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
    return frequency

class SequenceSet(Dataset):
    def __init__(self, path, k = 2):
        self.scores = []
        self.predictors = []
        with open(path) as f:
            for line in f:
                sequence_1, sequence_2, score = line.strip().split("\t")[:3]
                if score == 0:
                    continue
                score = -float(score)/10
                #predictor = list(count_frequency(sequence_1,k)) + [len(sequence_1)/100] + list(count_frequency(sequence_2,k)) + [len(sequence_2)/100]
                predictor = list(count_frequency(sequence_1,k)) + list(count_frequency(sequence_2,k)) 
                predictor = torch.tensor(predictor).double()
                self.scores.append(score)
                self.predictors.append(predictor)
        self.scores = torch.tensor(self.scores).double()
        
    def __len__(self):
        return min(len(self.scores),200000)
    
    def __getitem__(self,idx):
        return self.predictors[idx], self.scores[idx]
            
class GEV(nn.Module):
    def __init__(self,k):
        super().__init__()
        self.linear_1  = nn.Linear(2*4**k, 512, bias=True)
        self.linear_2  = nn.Linear(512, 256, bias=True)
        self.linear_3  = nn.Linear(256, 64, bias=True)
        self.linear_4  = nn.Linear(64, 32, bias=True)
        self.linear_5  = nn.Linear(32, 2, bias=True)
    def forward(self, x, y = None):
        x = self.linear_1(x).relu()
        x = self.linear_2(x).relu() 
        x = self.linear_3(x).relu()
        x = self.linear_4(x).relu()
        x = self.linear_5(x)
        loc= x[:,0]
        scale = x[:,1].relu().clamp(0.01)        
        if y is not None:
            error = (np.euler_gamma*scale + loc - y).abs()/y.clamp(0.1)
            y = (y-loc)/scale
            return (y+torch.exp(-y)+torch.log(scale)).mean(), error.mean()
        else:
            return loc, scale


def main():
    parser = argparse.ArgumentParser(description='train background model')
    parser.add_argument('--train', '-t', type=str, required=True, help='training dataset')
    parser.add_argument('--validation', '-v', type=str, required=True, help='validation dataset')
    parser.add_argument('--models', '-m', type=str, required=True, help='directory to save model')
    parser.add_argument('--performance', '-p', type=str, required=True, help='training log')
    parser.add_argument('--device', '-d', type=str, default="cuda:0", choices=["cuda:0","cuda:1"], help='device to use')
    parser.add_argument('--word-size', '-w', type=int, default=2, help='word size to use')
    parser.add_argument('--learning-rate', '-lr', type=int, default=0.00001, help='learning rate to use')
    args = parser.parse_args()    
    device = args.device
    logger.info("Load training data ...")
    train_dataset = SequenceSet(args.train,args.word_size)
    train_sampler = RandomSampler(train_dataset)
    train_loader = DataLoader(train_dataset, sampler=train_sampler, batch_size=512)    
    logger.info("Load validation data ...")
    val_dataset = SequenceSet(args.validation,args.word_size)
    val_sampler = RandomSampler(val_dataset)
    val_loader = DataLoader(val_dataset, sampler=val_sampler, batch_size=512)
    logger.info("Initialize model ...")
    model = GEV(args.word_size).to(device)
    model = model.double()
    optimizer = Adam(model.parameters(), lr=args.learning_rate, weight_decay=1e-5)
    if not os.path.exists(args.models):
        os.mkdir(args.models)

    logger.info("Start training ...")
    fp = open(args.performance,"w")
    for i in range(10000):
        losses = []
        errors = []
        model.train()
        for X, y in train_loader:
            X, y = X.to(device), y.to(device)
            optimizer.zero_grad()
            nll, error = model(X, y)
            nll.backward()
            optimizer.step()
            #print(nll.item(),mse.item())
            losses.append(nll.item())
            errors.append(error.item())
        print("train", i,np.mean(losses),np.mean(errors),sep="\t")
        print("train", i,np.mean(losses),np.mean(errors),sep="\t",file=fp)
        losses = []
        errors = []
        model.eval()
        for X, y in val_loader:
            X, y = X.to(device), y.to(device)
            nll, error = model(X, y)
            losses.append(nll.item())
            errors.append(error.item())
        print("val", i,np.mean(losses),np.mean(errors),sep="\t")
        print("val", i,np.mean(losses),np.mean(errors),sep="\t",file=fp)
        #print(i,np.mean(losses),np.mean(np.sqrt(mses)),sep="\t",file=fp)
        fp.flush()
        torch.save(model.state_dict(),f"{args.models}/gumbel.background.{i}.pt")        

    fp.close()
if __name__ == "__main__":
    main()
