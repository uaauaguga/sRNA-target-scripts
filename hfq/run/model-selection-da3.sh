#!/bin/bash
mkdir -p output/scores/model.selection.models.b10.c64.da3
for i in {21..30};do
  echo $i
  [ -s output/scores/model.selection.models.b10.c64.da3/${i}.bed ] || scripts/inference.py  -f dataset/inference/GTDB.proteobacteria.fa -o output/scores/model.selection.models.b10.c64.da3/${i}.bed  -m models.b10.c64.da3/${i}.pt -d cuda:0
  scripts/hfq-enrichment.py -i output/scores/model.selection.models.b10.c64.da3/${i}.bed -o output/scores/model.selection.models.b10.c64.da3/${i}.txt
done
