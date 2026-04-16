#!/bin/bash
for genome_id in $(cat  genome-ids.txt);do
  #[ -s genomes/leaders/${genome_id}.bed ] || 
  scripts/extract-leader-ivs.py -i genomes/bed/${genome_id}.bed -o genomes/leaders/${genome_id}.bed -l 200 -r 100
  #[ -s output/leader2peak/${genome_id}.bed ] || 
  bedtools intersect -loj -a genomes/leaders/${genome_id}.bed -b output/peaks.by.genome/${genome_id}.sorted.merged.bed -s > output/leader2peak/${genome_id}.bed 
done
