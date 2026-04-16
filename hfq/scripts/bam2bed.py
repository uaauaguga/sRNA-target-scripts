#!/usr/bin/env python
import argparse
import os
import subprocess
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('convert')

def main():
    parser = argparse.ArgumentParser(description='convert paired end bam file to merged bed file')
    parser.add_argument('--bam','-i',type=str, required=True, help='input bam file')
    parser.add_argument('--bed','-o',type=str, required=True, help="output bed file")
    parser.add_argument('--strand','-s',type=str, required=True, choices = ["forward","reverse","no"], help="strand")
    parser.add_argument('--min-quality','-mq',type=int, default=10, help="Minimal mapping quality required")
    parser.add_argument('--layout','-l',type=str, default="paired", help="Paired end or single end")
    args = parser.parse_args()
    
    logger.info("convert bam file to bedpe format ...")    
    ftmp = open(args.bed + ".tmp","w")
    if args.layout == "paired":
        cmd = ["bedtools","bamtobed","-i",args.bam,'-bedpe','-mate1']
        print(" ".join(cmd))
        subprocess.run(cmd,stdout=ftmp,stderr=subprocess.DEVNULL)
        ftmp.close()
    else:
        cmd = ["bedtools","bamtobed","-i",args.bam]
        print(" ".join(cmd))
        subprocess.run(cmd,stdout=ftmp,stderr=subprocess.DEVNULL)
    fout =  open(args.bed,"w")
    if args.layout == "paired":
        logger.info("merged paired read ...")
        n_cross_strand = 0
        n_chimeric = 0
        n_multiple_mapped = 0
        n_passed = 0
        with open(args.bed + ".tmp") as f:
            for line in f:
                seq_id_1, s1, e1, seq_id_2, s2, e2, name, value, strand_1, strand_2 = line.strip().split("\t") 
                s1, e1 = int(s1), int(e1)
                s2, e2 = int(s2), int(e2)
                s = min(s1,s2)
                e = max(e1,e2)
                if seq_id_1 != seq_id_2:
                    n_chimeric += 1
                    continue
                seq_id = seq_id_1
                if strand_1 == strand_2:
                    n_cross_strand += 1
                    continue
                if args.strand == "reverse":
                    strand = strand_2
                else:
                    strand = strand_1
                value = int(value)
                if value < args.min_quality:
                    n_multiple_mapped += 1
                    continue
                n_passed += 1
                print(seq_id, s, e , name, value, strand, sep="\t", file=fout)
        logger.info(f"{n_cross_strand} cross strand reads.")
        logger.info(f"{n_chimeric} chimeric reads.")
        logger.info(f"{n_multiple_mapped} multiple mapped reads.")
        logger.info(f"{n_passed} passed reads.")
    else:
        logger.info("extract reads ...")
        n_multiple_mapped = 0
        n_passed = 0
        with open(args.bed + ".tmp") as f:
            for line in f:
                seq_id, start, end, name, value, strand = line.strip().split("\t")[:6]
                start, end = int(start), int(end)
                value = int(value)
                if value < args.min_quality:
                    n_multiple_mapped += 1
                    continue
                n_passed += 1
                if args.strand == "reverse":
                    strand = {"+":"-","-":"+"}[strand]
                print(seq_id, start, end, name, value, strand, sep="\t", file=fout)
        logger.info(f"{n_multiple_mapped} multiple mapped reads.")
        logger.info(f"{n_passed} passed reads.")
    fout.close()
    logger.info("remove tmp file ...")
    os.remove(args.bed + ".tmp")
    logger.info("all done.")


if __name__ == "__main__":
    main()
