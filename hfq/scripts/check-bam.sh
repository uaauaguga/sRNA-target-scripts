#!/bin/bash
bam=$1
bed=$2
samtools view -h -f 66 ${bam} | bedtools bamtobed -i - | head -n 100000 | bedtools intersect -loj -a - -b ${bed} | \
         cut -f 6,10,12 | awk '$1!="."{{print $2,int($1==$3)}}' | sort | uniq -c | sort -k 1 -nr | head -n 1000 
