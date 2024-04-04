#!/usr/bin/env python
import logging
import argparse
import pickle
import torch
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] [%(levelname)s] %(message)s')
logger = logging.getLogger("convert weights to numpy array")

def main():
    parser = argparse.ArgumentParser(description='convert model')
    parser.add_argument('--input', '-i', type=str, required=True, help='input model weights')
    parser.add_argument('--output', '-o', type=str, required=True, help='output model weights')
    args = parser.parse_args()


    logger.info("load weights ...")
    weights = torch.load(args.input,map_location=torch.device('cpu'))
    
    for key in weights:
        weights[key] = weights[key].numpy()

    logger.info("saving converted weights ...")
    fout = open(args.output,"wb")
    pickle.dump(weights, fout)
    fout.close()

    logger.info("all done .")
    

if __name__ == "__main__":
    main()
