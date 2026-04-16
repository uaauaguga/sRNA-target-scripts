#!/usr/bin/env python
import argparse

def main():
    parser = argparse.ArgumentParser(description='pick confident sites')
    parser.add_argument('--input','-i',type=str, required=True, help="input scores in bed format")
    parser.add_argument('--output','-o',type=str, required=True, help="selected sites in bed format")
    parser.add_argument('--min-score','-ms',type=float,default=0.5,help="min score to consider")
    parser.add_argument('--min-distance','-md',type=int,default=50,help="two sites should have a distance >= this value")
    args = parser.parse_args()

    fout = open(args.output,"w")
    with open(args.input) as fin:
        local_max_score = -1
        local_max_position = -1
        last_chrom_id, last_start, last_end  = "", -1000, -1000
        for line in fin:
            fields = line.strip().split("\t")
            chrom_id, start, end, _, score, strand = fields[:6]
            start, end, score = int(start), int(end), float(score)
            # ignore positions with low score
            if score < args.min_score:
                continue
            # update current local statistics
            if chrom_id == last_chrom_id and start - last_end < args.min_distance:
                if score > local_max_score:
                    local_max_score = score
                    local_max_position = int((start+ end)/2)
                    local_chrom_id = chrom_id
            else:
                # save current local statistics, if any
                if local_max_position >= 0:
                    print(local_chrom_id,local_max_position,local_max_position+1,".",local_max_score,sep="\t",file=fout)
                # start a new local statistics
                local_max_position = int((start+end)/2)
                local_max_score = score
                local_chrom_id = chrom_id
            last_chrom_id, last_start, last_end = chrom_id, start, end


if __name__ == "__main__":
    main()
