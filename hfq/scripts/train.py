#!/usr/bin/env python
import os
import torch
from torch.optim import Adam
from dataset import SequenceSet,collate_fn
from torch.utils.data import DataLoader, SequentialSampler,RandomSampler
from torch.functional import F
import numpy as np
from sklearn.metrics import precision_recall_curve, auc, roc_auc_score, roc_curve
from model import CNNClassifier
import argparse
import pandas as pd
from scipy.interpolate import interp1d
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('train promoter predicter')
from ushuffle import shuffle

def main():
    parser = argparse.ArgumentParser(description='Train Sequence Classifier')
    parser.add_argument('--device','-d',default="cuda:1",choices=["cuda:0","cuda:1","cpu"],help="Device to run the model")
    parser.add_argument('--train-positive','-tp',required=True,help="positive sequence for training")
    parser.add_argument('--train-negative','-tn',required=True,help="negative sequence for training")
    parser.add_argument('--val-positive','-vp',required=True,help="positive sequence for evaluation")
    parser.add_argument('--val-negative','-vn',required=True,help="negative sequence for evaluation")
    parser.add_argument('--models','-m',required=True,help="directory to save models")
    parser.add_argument('--n-channels', '-c', type=int, default=64, help="number of channels to use")
    parser.add_argument('--n-blocks', '-b', type=int, default=5, help="number of resnet blocks to use")
    parser.add_argument('--kernel-size', '-k', type=int, default=5, help="convolution kernel size in the res-block")
    parser.add_argument('--length', '-l', type=int, default=100, help="length of the sequence set")
    parser.add_argument('--positive-fraction', '-pf', type=float, default=0.1, help="positive fraction")
    parser.add_argument('--shuffled-fraction', '-sf', type=float, default=0.0, help="shuffled sequence fraction in negative sequence")
    parser.add_argument('--batch-size', '-bs', type=int, default=512, help="batch size for ttraining and evaluation")
    parser.add_argument('--unstratify', '-us', action = "store_true",  help="whether stratify the sample when sampling")
    parser.add_argument('--domain-adversarial', '-da', type=float, default = 0.5, help="weight of domain adversarial loss")
    parser.add_argument('--crop-fraction', '-cf', type=float, default = 0.05, help="fraction to crop")
    args = parser.parse_args()

    if not os.path.exists(args.models):
        os.mkdir(args.models)

    train_set = SequenceSet(args.train_positive, args.train_negative, 
                           length=args.length, positive_fraction=args.positive_fraction, crop_fraction = args.crop_fraction,
                           shuffled_fraction = args.shuffled_fraction, stratified = not args.unstratify)
    if args.domain_adversarial > 0:
        n_domains = len(train_set.group_ids)
    else:
        n_domains = None
    train_sampler = RandomSampler(train_set)
    train_loader = DataLoader(train_set, sampler=train_sampler, batch_size=args.batch_size, collate_fn=collate_fn)

    val_set = SequenceSet(args.val_positive, args.val_negative, length=args.length, 
                          positive_fraction=args.positive_fraction, shuffled_fraction = 1,  stratified =  not args.unstratify)
    val_sampler = RandomSampler(val_set)
    val_loader = DataLoader(val_set, sampler=val_sampler, batch_size=args.batch_size, collate_fn=collate_fn)
    
    model = CNNClassifier(n_domains = n_domains, n_channels=args.n_channels, 
                          kernel_size = args.kernel_size, n_blocks = args.n_blocks)
    model = model.to(args.device)
    optimizer = Adam(model.parameters(), lr=0.0001)

    cls_criterion = torch.nn.CrossEntropyLoss()
    if args.domain_adversarial:
        da_criterion = torch.nn.CrossEntropyLoss()
    train_losses = []
    da_losses = []
    for e in range(1000):
        i = 0
        for sequences, sequence_labels, group_labels in train_loader:
            sequences, sequence_labels = sequences.to(args.device), sequence_labels.to(args.device)
            group_labels = group_labels.to(args.device)
            optimizer.zero_grad()
            if args.domain_adversarial > 0:
                cls_logits, da_logits = model(sequences)                
            else:
                cls_logits  = model(sequences)
            accuracy = (cls_logits.argmax(axis=1)==sequence_labels).sum()/sequence_labels.shape[0]
            cls_loss = cls_criterion(cls_logits,sequence_labels)
            if args.domain_adversarial > 0:
                # add domain adversarial loss
                da_loss = da_criterion(da_logits, group_labels)
                loss = args.domain_adversarial*args.domain_adversarial + cls_loss
            else:
                loss = cls_loss
            loss.backward()
            optimizer.step()
            train_losses.append(cls_loss.item())
            if args.domain_adversarial > 0:
                da_losses.append(da_loss.item())      
            if i%256 == 0:
                model = model.eval()
                if args.domain_adversarial:
                    logger.info(f"{e}-{i} train: {np.mean(train_losses)} {np.mean(da_losses)} {accuracy.item()}")
                else:
                    logger.info(f"{e}-{i} train: {np.mean(train_losses)} {accuracy.item()}")
                train_losses = []
                da_losses = []
                records = []
                correct, N = 0, 0
                losses = []
                j = 0
                y_trues, pred_probas = [], []
                for sequences, sequence_labels, group_labels in val_loader:
                    if j > 100:
                        break
                    sequences, sequence_labels = sequences.to(args.device), sequence_labels.to(args.device)
                    group_labels = group_labels.to(args.device)
                    if args.domain_adversarial > 0:
                        cls_logits, da_logits = model(sequences)
                    else:
                        cls_logits  = model(sequences)
                    cls_loss = cls_criterion(cls_logits,sequence_labels)
                    loss = cls_loss
                    losses.append(cls_loss.item())
                    pred_proba = torch.softmax(cls_logits,-1)[:,1].cpu().detach().numpy()
                    #best_index = np.argmax(pred_proba)
                    #if (pred_proba[best_index] > 0.75) and (sequence_labels[best_index] == 1):
                    #    s = ""
                    #    for loc in range(sequences[best_index].shape[1]):
                    #        c = torch.where(sequences[best_index,:,loc]==1)[0].item()
                    #        s += "ACGT"[c]
                    #    print(pred_proba[best_index], s)
                    y_true = sequence_labels.cpu().detach().numpy()
                    y_trues += list(y_true)
                    correct += ((pred_proba>0.5).astype(int)==y_true).sum()
                    N  += pred_proba.shape[0]
                    pred_probas += list(pred_proba)
                fpr, tpr, thresholds = roc_curve(y_trues, pred_probas)
                AUROC = auc(fpr, tpr)
                recall05 = round(float(interp1d(fpr,tpr)(0.05)),3)
                precision, recall, thresholds = precision_recall_curve(y_trues, pred_probas)
                AUPRC = auc(recall,precision)
                logger.info(f"{e}-{i} val: {np.mean(losses)} {correct/N} {AUPRC} {AUROC} {recall05}")
                model = model.train()
            i += 1
        logger.info(f"saving model at epoch {e} ...")
        torch.save(model.state_dict(),f"{args.models}/{e}.pt")
if __name__ == "__main__":
    main()
