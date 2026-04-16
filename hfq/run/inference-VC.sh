#!/bin/bash
pt=models.b10.c64.with.targets.da0.5/36.pt
clade=VC
indir=dataset/genomes/$clade
outdir=output/predictions/$clade
mkdir -p $outdir
for fasta in $(ls $indir);do
  genome_id=${fasta%.*}
  fasta=$indir/${genome_id}.fa
  bed=$outdir/${genome_id}.bed
  log=$outdir/${genome_id}.log
  [ -s $bed ] || scripts/inference.py -f $fasta -o $bed --n-blocks 10 -c 0 -rc -m $pt > $log 2>&1
done
