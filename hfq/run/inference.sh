#!/bin/bash
m=da2.only.target.36
pt=models.b10.c64.with.targets.da0.5/36.pt
#pt=models.b10.c64.only.targets.da2/36.pt
mkdir -p output/assemblies.${m}
for genome_id in $(cat entero9.txt);do
  fasta=assemblies/${genome_id}.fa
  bed=output/assemblies.${m}/${genome_id}.bed
  log=output/assemblies.${m}/${genome_id}.log
  [ -s $bed ] || scripts/inference.py -f $fasta -o $bed --n-blocks 10 -c 0 -rc -m $pt > $log 2>&1
done
