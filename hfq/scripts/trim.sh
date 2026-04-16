#!/bin/bash
dataset=$1
run_id=$2
if [ ! -s output/trimmed/$dataset/${run_id}_1.fastq.gz ];then
  trim_galore --cores 4 --length 10 --paired -a AGATCGGAAGAGCACACGTCTGAA -a2 AGATCGGAAGAGCGTCGTG --clip_R1 9 --three_prime_clip_R2 9 --output_dir output/trimmed/$dataset data/$dataset/${run_id}_1.fastq.gz data/$dataset/${run_id}_2.fastq.gz > output/trimmed/$dataset/${run_id}.log 2>&1
  mv output/trimmed/$dataset/${run_id}_1_val_1.fq.gz output/trimmed/$dataset/${run_id}_1.fastq.gz
  mv output/trimmed/$dataset/${run_id}_2_val_2.fq.gz output/trimmed/$dataset/${run_id}_2.fastq.gz
fi
