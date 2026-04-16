from collections import defaultdict
import numpy as np
import os
import torch
from torch.utils.data import Dataset
import logging
import random
from tqdm import tqdm
import sys
import re
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] %(message)s')
logger = logging.getLogger("Load training instances")
from ushuffle import shuffle


def collate_fn(examples):
    sequences = []
    sequence_labels = []
    group_labels = []
    for sequence, sequence_label, group_label in examples:
        sequences.append(sequence)
        sequence_labels.append(sequence_label)
        group_labels.append(group_label)
    sequence_labels = torch.tensor(sequence_labels)
    group_labels = torch.tensor(group_labels)
    # batch , channel , length
    sequences = torch.cat(sequences)
    return sequences, sequence_labels, group_labels

DNA_alp = dict(zip("ACGT",list(range(4))))
protein_alp = dict(zip("ACDEFGHIKLMNPQRSTVWY",list(range(20))))

def onehot(sequence,alp="DNA"):
    assert alp in ["DNA","protein"], "the alphabet should be either DNA or protein"
    if alp == "DNA":
        c2i = DNA_alp
        x = torch.zeros((4,len(sequence)))
    else:
        c2i = protein_alp 
        x = torch.zeros((20,len(sequence)))
    tokens = []
    indices = []
    for i,c in enumerate(sequence):
        if c in c2i:
            tokens.append(c2i[c])  
            indices.append(i)
    x[tokens,indices] = 1
    return x[None,:,:]

def unify(sequence, L):
    if len(sequence) < L:
        sequence = sequence + "N"*(L-len(sequence))
    elif len(sequence) > L:
        hL = int(L/2)
        center = int(len(sequence)/2)
        sequence = sequence[center-hL:center+hL]
    return sequence                

class SequenceSet(Dataset):
    """
    load  sequences
    input:
       positive sequence
       background sequence
    """

    def __init__(self,sequences,background, positive_fraction=0.5, 
                 shuffled_fraction = 0, length=50, stratified = False, crop_fraction = 0):
        self.positive_fraction = positive_fraction
        self.shuffled_fraction = shuffled_fraction
        self.crop_fraction = crop_fraction
        self.length = length
        logger.info("Load positive sequences ...")
        self.positives = defaultdict(list)
        with open(sequences) as f:
            for header in f:
                sequence = next(f).strip().upper().replace("U","T")
                #if len(sequence) != length:
                #    continue
                sequence = unify(sequence, length)
                seq_id = header[1:].strip().split(" ")[0]
                if stratified:
                    group_id = seq_id.split(":")[0]             
                else:
                    group_id = "dummy"
                self.positives[group_id].append((seq_id, sequence))
        self.background = defaultdict(list)
        logger.info("Load background sequences ...")
        with open(background) as f:
           for header in f:
                sequence = next(f).strip().upper().replace("U","T")
                #if len(sequence) != length:
                #    continue
                sequence = unify(sequence, length)
                seq_id = header[1:].strip().split(" ")[0]
                if stratified:
                    group_id = seq_id.split(":")[0]
                else:
                    group_id = "dummy"
                self.background[group_id].append((seq_id, sequence))
        self.group_ids = [ group_id for group_id in self.background]
        self.common_group_ids = [ group_id for group_id in self.background if group_id in self.positives] 
        self.g2i = {}
        for i, group_id in enumerate(self.group_ids):
            self.g2i[group_id] = i
        logger.info(f"{len(self.group_ids)} groups in background set .")
        logger.info(f"{len(self.common_group_ids)} present in both positive set and background set .")
        

    def __len__(self):
        return 10000

    def __getitem__(self,idx):
        group_id = self.group_ids[idx%len(self.group_ids)]
        if (random.random() < self.positive_fraction) and (len(self.positives[group_id]) > 0):
            seq_id, sequence = self.positives[group_id][np.random.randint(0,len(self.positives[group_id]))]
            sequence_label = 1
        else:
            seq_id, sequence = self.background[group_id][np.random.randint(0,len(self.background[group_id]))]
            if random.random() < self.shuffled_fraction:
                sequence = shuffle(sequence,len(sequence),2)
            sequence_label = 0
        group_id = seq_id[:seq_id.find(":")]
        sequence = onehot(sequence)
        if np.random.rand() < self.crop_fraction:
            cropped_location = np.random.randint(int(self.length/2))
            sequence[:,:,cropped_location:] = 0
        group_label = self.g2i[group_id]
        return sequence, sequence_label, group_label
