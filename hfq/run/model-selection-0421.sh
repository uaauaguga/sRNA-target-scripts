#!/bin/bash
params=2025.04.models.b10.c64.with.targets.da0.5
mkdir -p output/scores/$params
for i in {20..30};do
  echo $i
  [ -s output/scores/$params/${i}.bed ] || scripts/inference.py  -f dataset/inference/GTDB.proteobacteria.fa -o output/scores/$params/${i}.bed -m $params/${i}.pt -d cuda:0
  scripts/hfq-enrichment.py -i output/scores/$params/${i}.bed -o output/scores/$params/${i}.txt
done
