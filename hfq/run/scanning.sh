#!/bin/bash
pt=models.b10.c64.with.targets.da0.5/36.pt
for fasta in $(ls dataset/inference/GTDB | grep '.fa$');do
  chunk_id=${fasta%.*}
  [ -s output/scores/GTDB/${chunk_id}.bed ] || scripts/inference.py  -f dataset/inference/GTDB/${chunk_id}.fa -o output/scores/GTDB/${chunk_id}.bed  -m  $pt > output/scores/GTDB/${chunk_id}.log 2>&1
  [ -s output/scores/GTDB/${chunk_id}.txt ] || scripts/hfq-enrichment.py -i output/scores/GTDB/${chunk_id}.bed -o output/scores/GTDB/${chunk_id}.txt >  output/scores/GTDB/${chunk_id}.auroc.log 2>&1
  scripts/hfq-Z-score.py -i output/scores/GTDB/${chunk_id}.bed -o output/scores/GTDB/${chunk_id}.Z-score.txt
done
