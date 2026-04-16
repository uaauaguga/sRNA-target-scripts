#!/usr/bin/env python
import argparse
from pyfaidx import Fasta

def main():
    parser = argparse.ArgumentParser(description='extract leader sequences')
    parser.add_argument('--input',  '-i', required=True, help="input CDS bed")
    parser.add_argument('--output', '-o', required=True, help="input leader fasta")
    parser.add_argument('--left',   '-l', type=int, default = 200,  help="left flanking length")
    parser.add_argument('--right',  '-r', type=int,default = 100,  help="right flanking length")    
    args = parser.parse_args()

    fout = open(args.output,"w")
    with open(args.input) as fin:
        for line in fin:
            seq_id, start, end, name, score, strand = line.strip().split("\t")[:6]
            start, end = int(start), int(end)
            if strand == "+":
                s = max(0,start - args.left)
                e = start + args.right
            else:
                s = max(0,end - args.right)
                e = end + args.left
            print(seq_id, start, end, name, ".", strand, sep="\t",file=fout)
    fout.close()


if __name__ == "__main__":
    main()
