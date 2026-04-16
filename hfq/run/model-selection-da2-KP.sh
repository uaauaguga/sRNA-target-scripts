#!/bin/bash
outdir=output/model-selection/KP/da1
mkdir -p $outdir
for i in {0..40};do
  echo $i
  [ -s $outdir/${i}.bed ] || scripts/inference.py  -f dataset/KP.test.fa -o $outdir/${i}.bed -m models.b10.c64.da2/${i}.pt 
  tail test.bed
  scripts/performance-evaluation.py -i $outdir/${i}.bed -o test.txt 
done
