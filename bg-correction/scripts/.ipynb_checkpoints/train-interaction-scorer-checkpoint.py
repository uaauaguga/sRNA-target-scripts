#!/usr/bin/env python
import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader, SequentialSampler,RandomSampler
from torch.functional import F
from torch import nn
from collections import defaultdict
from torch.optim import Adam
import sys
import argparse
import os
from torch_scatter import scatter
from model import InteractionScorer
from tqdm import tqdm
from sklearn.metrics import roc_auc_score,roc_curve
import logging
from scipy.interpolate import interp1d
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] %(message)s')
logger = logging.getLogger("model training")

from itertools import product

def get_rc(s):
    comp_lut = {"A":["T"],"C":["G"],"G":["C","T"],"T":["A","G"]}
    last_comps = [""]
    for c in s:        
        comps = []
        for comp in last_comps:
            for cc in comp_lut[c]:
                comps.append(comp + cc)
        last_comps = comps
    return [comp[::-1] for comp in comps]   

k = 6
kmers = []
for i,kmer in enumerate(list(product(*["ACGT"]*k))):
    kmer = "".join(kmer)
    kmers.append(kmer)

kmer_rcs = {}
for kmer in kmers:
    kmer_rcs[kmer] = get_rc(kmer)

DNA_alp = dict(zip("ACGT",list(range(4))))
def onehot(sequence,length):
    c2i = DNA_alp
    x = torch.zeros((4,length))
    tokens = []
    indices = []
    for i,c in enumerate(sequence[:length]):
        if c in c2i:
            tokens.append(c2i[c])
            indices.append(i)
    x[tokens,indices] = 1
    return x.float()

kmer_hybrid_energies = {}
kmer_hybrid_energy_cutoff = -0.1

def select_candidate(s1,s2,flanking=40,min_distance=32):
    kmer_rc_set = {}
    for i in range(len(s1)-k):
        kmer = s1[i:i+k]
        if kmer not in kmer_rcs:
            continue
        for kmer_rc in kmer_rcs[kmer]:
            if kmer_rc not in kmer_rc_set:
                kmer_rc_set[kmer_rc] = []
            kmer_rc_set[kmer_rc].append(i)
    n = 0
    last_j = -100
    last_energy = 0
    flanked_s1 = "N"*flanking + s1 + "N"*flanking
    flanked_s2 = "N"*flanking + s2 + "N"*flanking    
    candidate_positions = []
    for j in range(len(s2)-k):
        kmer = s2[j:j+k]
        if kmer in kmer_rc_set: 
            # reverse complementary of current kmer present in first sequence
            for i in kmer_rc_set[kmer]:
                # extract correspnding positions in first sequence
                energy = kmer_hybrid_energies[(s1[i:i+k], kmer)]
                if  energy < kmer_hybrid_energy_cutoff:
                    if energy < last_energy:
                        last_energy = energy
                        candidate_position = (i,j)
                    if j - last_j > min_distance:
                        if candidate_position is not None:
                            candidate_positions.append(candidate_position)
                        last_energy = 0
                        last_j = j
                        candidate_position = None
    candidates_1, candidates_2 = [], []
    for i,j in candidate_positions:
        candidate_1 = flanked_s1[i:i+k+2*flanking]
        candidate_2 = flanked_s2[j:j+k+2*flanking]
        candidate_1 = onehot(candidate_1,k+2*flanking)                    
        candidate_2 = onehot(candidate_2,k+2*flanking)
        candidates_1.append(candidate_1)
        candidates_2.append(candidate_2)
    L = 2*flanking+k  
    if len(candidates_1) > 0:
        candidates_1 = torch.stack(candidates_1)[:,:,:,None].repeat([1,1,1,L])
        candidates_2 = torch.stack(candidates_2)[:,:,None,:].repeat([1,1,L,1])
        outer = torch.cat([candidates_1,candidates_2],axis=1)
    else:
        outer = torch.zeros(1,8,L,L)
    return outer

class SequencePairSet(Dataset):
    def __init__(self, path, max_length=512,chunk_size=5000):
        self.labels = []
        self.pairs = []
        self.handle = open(path)
        self.chunk_size = chunk_size
        self.max_length = max_length
        self.update()
              
    def __len__(self):
        return self.chunk_size
    
    def __getitem__(self,idx):        
        return self.pairs[idx], self.probs[idx]
    
    def update(self):
        logger.info("Update dataset ...")
        i = 0
        self.pairs = []
        self.probs = []
        n_candidates = 0
        while True:
            line = self.handle.readline()
            if not line:
                self.handle.seek(0)
                line = self.handle.readline()            
            sequence_1, sequence_2, energy, prob = line.strip().split("\t")[:4]
            #prob = float(prob)
            energy = float(energy)
            #prob = 1 - np.exp(float(energy)/13)
            prob = int(energy < -9.20188)
            if (len(sequence_1) > self.max_length) or (len(sequence_2) > self.max_length):
                continue
            i += 1
            self.probs.append(prob)
            candidates = select_candidate(sequence_1,sequence_2,flanking=37,min_distance=40)
            n_candidates += len(candidates)
            self.pairs.append(candidates)      
            if i == self.chunk_size:
                break
        logger.info(f"{int(n_candidates/len(self.pairs))} candidates per-sequence are selected.")

                
class RNAPairSet(Dataset):
    def __init__(self, path, cutoff = 1, hard_fraction=0.5,
                 max_length=512, chunk_size=5000, positive_fraction=0.2,clip=False):
        self.max_length = max_length
        self.hard_fraction = hard_fraction
        self.positive_fraction = positive_fraction
        self.clip = clip
        self.sRNAs = {}
        self.targets = {}        
        self.hard_positive_pairs = [] 
        self.hard_positive_pairs_tip_scores = []
        self.easy_positive_pairs = [] 
        self.easy_positive_pairs_tip_scores = []
        self.hard_negative_pairs = [] 
        self.hard_negative_pairs_tip_scores = []
        self.easy_negative_pairs = [] 
        self.easy_negative_pairs_tip_scores = []
        group_idx = 0
        for prefix in open(path).read().strip().split("\n"):
            if (not os.path.exists(f"{prefix}.1.fa")) or (not os.path.exists(f"{prefix}.2.fa")) or (not os.path.exists(f"{prefix}.txt")):
                continue
            logger.info(f"Load data from {prefix} ...")
            group = str(group_idx).zfill(4)
            logger.info(group)
            with open(f"{prefix}.1.fa") as f:
                for header in f:
                    seq_id = header[1:].strip().split(" ")[0]
                    sequence = next(f).strip()
                    if len(sequence) > self.max_length:
                        continue
                    self.sRNAs[group + ":" + seq_id] = sequence
            logger.info("Load targets ...")
            with open(f"{prefix}.2.fa") as f:
                for header in f:
                    seq_id = header[1:].strip().split(" ")[0]
                    sequence = next(f).strip()
                    if len(sequence) > self.max_length:
                        continue
                    self.targets[group + ":" + seq_id] = sequence              
            logger.info("Load interactions ...")
            with open(f"{prefix}.txt") as f:
                for line in f:
                    fields = line.strip().split("\t")
                    sRNA_id, target_id, tip_score, short_score, long_score = fields
                    tip_score, min_score, max_score = float(tip_score), float(short_score), float(long_score)
                    if (group + ":" + sRNA_id) not in self.sRNAs:
                        continue
                    if (group + ":" + target_id) not in self.targets:
                        continue                    
                    if max_score < cutoff:
                        label = 0
                    elif min_score > cutoff:
                        label = 1
                    else:
                        # not confident
                        continue
                    tip_label = int(tip_score > cutoff)
                    if label == tip_label:
                        # easy instance
                        if label == 1:                        
                            self.easy_positive_pairs.append((group + ":" + sRNA_id, group + ":" + target_id))
                            self.easy_positive_pairs_tip_scores.append(tip_score)
                        else:
                            self.easy_negative_pairs.append((group + ":" + sRNA_id, group + ":" + target_id))
                            self.easy_negative_pairs_tip_scores.append(tip_score)
                    else:
                        # hard instance
                        if label == 1:
                            self.hard_positive_pairs.append((group + ":" + sRNA_id, group + ":" + target_id))
                            self.hard_positive_pairs_tip_scores.append(tip_score)
                        else:
                            self.hard_negative_pairs.append((group + ":" + sRNA_id, group + ":" + target_id))
                            self.hard_negative_pairs_tip_scores.append(tip_score)                                                    
            group_idx += 1
        logger.info(f"{len(self.easy_positive_pairs)} easy positive pairs, ")
        logger.info(f"{len(self.hard_positive_pairs)} hard positive pairs, ")
        logger.info(f"{len(self.easy_negative_pairs)} easy negative pairs,")
        logger.info(f"{len(self.hard_negative_pairs)} hard negative pairs loaded.")
        logger.info(f"{len(self.sRNAs)} {len(self.targets)}")
        self.chunk_size = chunk_size
        self.max_length = max_length
        self.update()
              
    def __len__(self):
        return self.chunk_size
    
    def __getitem__(self,idx):
        idx = idx%len(self.pairs)
        return self.pairs[idx], self.labels[idx]
    
    def update(self):
        logger.info("Update dataset ...")
        self.labels = []
        self.pairs = []
        self.scores = []
        n_candidates = 0
        i = 0
        counter = defaultdict(int)
        #for i in range(self.chunk_size):
        j = 0
        while i < self.chunk_size:
            j += 1
            # determine the label
            label = int(np.random.rand() < self.positive_fraction)
            if np.random.rand() < self.hard_fraction:
                # sampling hard instance
                if label == 1:
                    idx = np.random.randint(len(self.hard_positive_pairs))
                    score = self.hard_positive_pairs_tip_scores[idx]
                    sRNA_id, target_id = self.hard_positive_pairs[idx]                 
                else:
                    idx = np.random.randint(len(self.hard_negative_pairs))
                    score = self.hard_negative_pairs_tip_scores[idx]
                    sRNA_id, target_id = self.hard_negative_pairs[idx]
            else:
                # sampling easy instance
                if label == 1:
                    idx = np.random.randint(len(self.easy_positive_pairs))
                    score = self.easy_positive_pairs_tip_scores[idx]
                    sRNA_id, target_id = self.easy_positive_pairs[idx]                 
                else:
                    idx = np.random.randint(len(self.easy_negative_pairs))
                    score = self.easy_negative_pairs_tip_scores[idx]
                    sRNA_id, target_id = self.easy_negative_pairs[idx]              
            group_id = sRNA_id.split(":")[0]
            if self.clip and (counter[group_id] > 500):
                continue
            counter[group_id] += 1
            sequence_1 = self.sRNAs[sRNA_id]
            sequence_2 = self.targets[target_id]
            self.labels.append(label)
            self.scores.append(score)
            candidates = select_candidate(sequence_1,sequence_2,flanking=37,min_distance=40)
            n_candidates += len(candidates)
            self.pairs.append(candidates) 
            i += 1
        if self.clip:    
            logger.info(f"{int(100*(j-i)/j)} of sampled instances are rejected .")
        AUROC = roc_auc_score(self.labels,self.scores)
        fpr, tpr, thresholds = roc_curve(self.labels,self.scores)  
        fpr2tpr = interp1d(fpr,tpr)
        recall01 = round(float(fpr2tpr(0.01)),3)
        recall05 = round(float(fpr2tpr(0.05)),3)
        recall10 = round(float(fpr2tpr(0.1)),3)
        logger.info(f"Baseline AUROC is {AUROC} .")
        logger.info(f"recalls: {recall01} {recall05} {recall10}")
        logger.info(f"{int(n_candidates/self.chunk_size)} candidates per-sequence are selected.")
                
                
                
def collate(bags):
    instances = []
    indices = []
    labels = []
    for i,(bag,label) in enumerate(bags):
        labels.append(label)
        for instance in bag:
            instances.append(instance[None,:,:,:])
            indices.append(i)
    instances = torch.vstack(instances)
    indices = torch.tensor(indices)
    #label by bags
    labels = torch.tensor(labels)
    return instances, indices, labels.float()

epsilon = 0.0000001
            
def main():
    parser = argparse.ArgumentParser(description='train background model')
    parser.add_argument('--train', '-t', type=str, required=True, help='training dataset')
    parser.add_argument('--validation', '-v', type=str, required=True, help='validation dataset')
    parser.add_argument('--kmer-hybrid-energies', '-khe', type=str, default="hexamer-energy.txt", help='hybrid energies of hexamers')
    parser.add_argument('--max-energy', '-me', type=float, default=-0.1, help='maximal of hybrid energy')
    parser.add_argument('--models', '-m', type=str, required=True, help='directory to save model')
    parser.add_argument('--performance', '-p', type=str, required=True, help='training log')
    #parser.add_argument('--pooling', type=str, choices=["max","noisy-or"], default="noisy-or",
    #                    help='method for pooling instance at bag level')
    parser.add_argument('--device', '-d', type=str, default="cuda:0", choices=["cuda:0","cuda:1"], help='device to use')
    parser.add_argument('--learning-rate', '-lr', type=float, default=0.0001, help='learning rate to use')
    parser.add_argument('--data', type=str, choices=["1","2"], help='data used for training')
    parser.add_argument('--resume-from', '-rf', type=str, help='resume model weights for training')
    parser.add_argument('--freeze', action="store_true", help='whether fix some model parameter')
    args = parser.parse_args()    
    global kmer_hybrid_energies
    global kmer_hybrid_energy_cutoff
    with open(args.kmer_hybrid_energies) as f:
        for line in f:
            hex1, hex2, energy = line.strip().split("\t")
            energy = float(energy)
            kmer_hybrid_energies[(hex1, hex2)] = energy
            kmer_hybrid_energies[(hex2, hex1)] = energy
    kmer_hybrid_energy_cutoff = args.max_energy
    device = args.device
    logger.info(f"Load training dataset from {args.train} ...")
    if args.data == "1":
        train_dataset = SequencePairSet(args.train)
    else:
        train_dataset = RNAPairSet(args.train,hard_fraction=0.5,positive_fraction=0.2,clip=True)
    train_sampler = RandomSampler(train_dataset)
    train_loader = DataLoader(train_dataset, sampler=train_sampler, batch_size=64,collate_fn=collate)  
    logger.info(f"Load validation dataset from {args.validation} ...")
    if args.data == "1":
        val_dataset = SequencePairSet(args.validation)
    else:
        val_dataset = RNAPairSet(args.validation, hard_fraction=0.5,positive_fraction=0.01)
    val_sampler = RandomSampler(val_dataset)
    val_loader = DataLoader(val_dataset, sampler=val_sampler, batch_size=64,collate_fn=collate)     
    logger.info("Intialize model ...")
    model = InteractionScorer()
    if args.resume_from is not None:
        logger.info(f"Initialize model parameters from {args.resume_from} ...")
        state_dict = torch.load(args.resume_from,map_location = args.device)
        model.load_state_dict(state_dict)
    model = model.to(device)
    if args.freeze:
        logger.info("Freeze part of model parameter ...")
        for n, p in model.named_parameters():
            #if n.split(".")[0] not in ["fc"]: #and (not n.startswith("res3.downsample")):
            if n.split(".")[0] not in ["fc","res3"]:
                p.requires_grad = False
    optimizer = Adam(model.parameters(), lr=args.learning_rate)
    if not os.path.exists(args.models):
        os.mkdir(args.models)
    logger.info("Start training ...")
    fp = open(args.performance,"w")
    for i in range(10000):
        train_losses = []
        logger.info(f"training epoch {i+1} ...")
        model.train()
        torch.save(model.state_dict(),f"{args.models}/model.{i}.pt")   
        y_true, y_pred = [], []
        for instances, indices, labels in train_loader:
            #print(labels[:10])
            instances, indices, labels = instances.to(device), indices.to(device), labels.to(device)
            optimizer.zero_grad()
            scores = model(instances)
            scores, indices = scores.view(-1), indices.view(-1)
            pooled_scores = scatter(scores, indices,  reduce="max")
            #if args.pooling == "max":
            #    pooled_scores = scatter(scores, indices,  reduce="max")
            #else: 
            #    pooled_scores = torch.ones(labels.shape[0]).to(device)
            #    pooled_scores = 1 - torch.exp(scatter(torch.log(1 - scores).clamp(epsilon), indices,  reduce="sum"))
            #    pooled_scores = pooled_scores.clamp(epsilon)
            loss = F.binary_cross_entropy(pooled_scores,labels)
            loss.backward()
            optimizer.step()
            train_losses.append(loss.item())
            y_true += list((labels.detach().cpu().numpy()>0.5).astype(int))
            y_pred += list(pooled_scores.detach().cpu().numpy())
        train_AUROC = roc_auc_score(np.array(y_true).astype(int),np.array(y_pred))
        model.eval()
        y_true, y_pred = [], []
        val_losses = []
        for instances, indices, labels in val_loader:
            instances, indices, labels = instances.to(device), indices.to(device), labels.to(device)
            scores = model(instances)
            scores, indices = scores.view(-1), indices.view(-1)
            pooled_scores = scatter(scores, indices,  reduce="max")
            #if args.pooling == "max":
            #    pooled_scores = scatter(scores, indices,  reduce="max")
            #else:
            #    pooled_scores = 1 - torch.exp(scatter(torch.log(1 - scores).clamp(epsilon), indices,  reduce="sum"))
            #    pooled_scores = pooled_scores.clamp(epsilon)
            loss = F.binary_cross_entropy(pooled_scores,labels)
            val_losses.append(loss.item())
            y_true += list((labels.detach().cpu().numpy()>0.5).astype(int))
            y_pred += list(pooled_scores.detach().cpu().numpy())
        val_AUROC = roc_auc_score(np.array(y_true).astype(int),np.array(y_pred))
        fpr, tpr, thresholds = roc_curve(np.array(y_true).astype(int),np.array(y_pred))  
        fpr2tpr = interp1d(fpr,tpr)
        recall01 = round(float(fpr2tpr(0.01)),3)
        recall05 = round(float(fpr2tpr(0.05)),3)
        recall10 = round(float(fpr2tpr(0.1)),3)        
        print(i,np.mean(train_losses),train_AUROC,np.mean(val_losses),val_AUROC,recall01,recall05,recall10,sep="\t")
        print(i,np.mean(train_losses),train_AUROC,np.mean(val_losses),val_AUROC,recall01,recall05,recall10,sep="\t",file=fp)
        train_dataset.update()
        #val_dataset.update()
        fp.flush() 
    fp.close()
if __name__ == "__main__":
    main()
