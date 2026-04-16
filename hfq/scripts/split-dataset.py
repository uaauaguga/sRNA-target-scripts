#!/usr/bin/env python
import argparse
from collections import defaultdict
import os
import logging
import numpy as np
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] [%(levelname)s] %(message)s')
logger = logging.getLogger("groupping sequences")

def main():
    parser = argparse.ArgumentParser(description='split dataset')
    parser.add_argument('--input', '-i', type=str, required=True, help='input sequences')
    parser.add_argument('--train','-t', type=str ,required=True, help="training data")
    parser.add_argument('--validation','-v',type=str,required=True,help="validation data")
    parser.add_argument('--train-fraction','-tf',type=float,default=0.9,help="fraction of training data")
    args = parser.parse_args()
    
    ft = open(args.train,"w")
    fv = open(args.validation,"w")
    with open(args.input) as f:
        for header in f:
            sequence = next(f)
            fout = ft if  np.random.rand() < args.train_fraction else fv
            fout.write(header)
            fout.write(sequence)
    ft.close()
    fv.close()


if __name__ == "__main__":
    main()

