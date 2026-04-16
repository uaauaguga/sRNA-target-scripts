#!/bin/bash
mkdir -p output/binding.sites.by.genome
for genome_id in $(cat genome-ids.txt) ;do
  [ -s output/binding-sites-by-genome/${dataset}.fa ] || scripts/extract-binding-site-sequences.py -i output/peaks.by.genome/${genome_id}.sorted.merged.bed -g genomes/fasta/${genome_id}.fa -o output/binding.sites.by.genome/${genome_id}.fa -p 0.0005
done
