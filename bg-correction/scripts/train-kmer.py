#!/usr/bin/env python
from itertools import product
import torch
from torch.utils.data import Dataset, DataLoader, SequentialSampler,RandomSampler
from torch.functional import F
from torch import nn
from collections import defaultdict
from torch.optim import Adam, lr_scheduler
import numpy as np
import sys
import argparse
import os
import re
from sklearn.metrics import roc_auc_score
import words
words.init(5)

def count_frequency(sequence,k=5):
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
    def __init__(self, path, k = 5, n_samples = None):
        self.scores = []
        self.predictors = []
        self.n_samples = n_samples
        with open(path) as f:
            j = 0
            for line in f:
                j += 1
                if j%5000 == 0:
                    print(f"{round(j/1000,2)} K entries loaded .")
                sequence_1, sequence_2, score = line.strip().split("\t")[:3]
                score = -float(score)/10
                predictor = list(count_frequency(sequence_1,k)) + list(count_frequency(sequence_2,k))
                predictor = torch.tensor(predictor).float()
                self.scores.append(score)
                self.predictors.append(predictor)
        self.scores = torch.tensor(self.scores).float()
    def __len__(self):
        return 100000 if self.n_samples is None else self.n_samples
    
    def __getitem__(self,idx):
        idx = idx%len(self.scores)
        return self.predictors[idx], self.scores[idx]
            
class GEV(nn.Module):
    def __init__(self,k):
        super().__init__()
        self.linear_1  = nn.Linear(2*4**k, 512, bias=True)        
        self.linear_2  = nn.Linear(512, 256, bias=True)
        self.linear_3  = nn.Linear(256, 128, bias=True)
        self.gumbel_1  = nn.Linear(128, 64, bias=True)
        self.gumbel_2  = nn.Linear(64, 2, bias=True)
        
    def forward(self, x):
        out = self.linear_3(self.linear_2(self.linear_1(x).relu()).relu()).relu()
        gumbel_params = self.gumbel_2(self.gumbel_1(out).relu())
        locs = gumbel_params[:,0]
        scales = gumbel_params[:,1].relu().clamp(0.1)        
        return locs, scales


def gumbel_loss(locs, scales, y):
    mse0 = (np.euler_gamma*scales + locs - y)**2
    mse = (locs - y)**2
    y = (y-locs)/scales
    return (y+torch.exp(-y)+scales).mean(), mse0.mean(), mse.mean()

def main():
    parser = argparse.ArgumentParser(description='train background model')
    parser.add_argument('--train', '-t', type=str, required=True, help='input training dataset')
    parser.add_argument('--validation', '-v', type=str, required=True, help='input validation dataset')
    parser.add_argument('--models', '-m', type=str, required=True, help='directory to save model')
    parser.add_argument('--performance', '-p', type=str, required=True, help='training log')
    parser.add_argument('--device', '-d', type=str, default="cuda:0", choices=["cuda:0","cuda:1"], help='device to use')
    parser.add_argument('--word-size', '-w', type=int, default=5, help='word size to use')
    parser.add_argument('--learning-rate', '-lr', type=int, default=0.001, help='learning rate to use')
    args = parser.parse_args()    
    device = args.device
    train_dataset = SequenceSet(args.train, args.word_size,2000000)
    train_sampler = RandomSampler(train_dataset)
    train_loader = DataLoader(train_dataset, sampler=train_sampler, batch_size=2048)    

    val_dataset = SequenceSet(args.validation, args.word_size,2000)
    val_sampler = RandomSampler(val_dataset)
    val_loader = DataLoader(val_dataset, sampler=val_sampler, batch_size=256)    

    model = GEV(args.word_size).to(device)
    optimizer = Adam(model.parameters(), lr=args.learning_rate, weight_decay=1e-4)
    scheduler = lr_scheduler.ExponentialLR(optimizer,gamma=0.99)
    if not os.path.exists(args.models):
        os.mkdir(args.models)
    fp = open(args.performance,"w")
    for i in range(500):
        losses = []
        mses = []
        mse0s = []
        model.train()
        for X, y in train_loader:            
            X, y = X.to(device), y.to(device)
            mask = X[:,:int(X.shape[1]/2)] > 0
            optimizer.zero_grad()
            locs, scales = model(X)
            energy_nll, mse0, mse = gumbel_loss(locs, scales, y)
            loss = energy_nll
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            losses.append(loss.item())
            mses.append(mse.item())
            mse0s.append(mse0.item())
        scheduler.step()
        print("train", i,np.mean(losses),np.mean(np.sqrt(mse0s)),np.mean(np.sqrt(mses)),sep="\t",file=fp)
        print("train", i,np.mean(losses),np.mean(np.sqrt(mse0s)),np.mean(np.sqrt(mses)),sep="\t")
        #print(i,np.mean(losses),np.mean(mses),sep="\t")
        energy_losses = []
        mses = []
        mse0s = []
        AUROCs = []
        model.eval()
        for X, y in val_loader:
            X, y = X.to(device), y.to(device)
            locs, scales = model(X)
            energy_nll, mse0, mse = gumbel_loss(locs, scales, y)
            y_pred = (np.euler_gamma*scales + locs).detach().cpu().numpy()
            y_true = (y.cpu().numpy() > 1).astype(int)
            AUROC = roc_auc_score(y_true, y_pred) 
            energy_losses.append(energy_nll.item())
            mses.append(mse.item())
            mse0s.append(mse0.item())
            AUROCs.append(AUROC)
        print("val",i,np.mean(energy_losses),np.mean(np.sqrt(mse0s)),np.mean(np.sqrt(mses)),np.mean(AUROCs),sep="\t",file=fp)
        print("val",i,np.mean(energy_losses),np.mean(np.sqrt(mse0s)),np.mean(np.sqrt(mses)),np.mean(AUROCs),sep="\t")
        fp.flush()
        if i%500 == 1: 
            torch.save(model.state_dict(),f"{args.models}/gumbel.background.{i}.pt")        
    fp.close()
if __name__ == "__main__":
    main()
