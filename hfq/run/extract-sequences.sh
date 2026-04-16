#!/bin/bash
while read dataset asm_id;do
  echo $dataset $asm_id
  [ -s output/binding-sites/${dataset}.fa ] || scripts/extract-binding-site-sequences.py -i output/peaks/${dataset}.bed -g genomes/fasta/${asm_id}.fa -o output/binding-sites/${dataset}.fa -p 0.0005
done < dataset2genome.txt #clip-dataset-genome.txt	
