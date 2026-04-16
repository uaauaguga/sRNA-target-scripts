#!/bin/bash
mkdir -p output/peaks.by.genome.shuffled
mkdir -p output/peaks.by.genome.shuffled.annotated
for genome_id in $(cat  genome-ids.txt);do
  bedtools shuffle -i output/peaks.by.genome/${genome_id}.sorted.merged.bed -g genomes/fasta/${genome_id}.fa.fai > output/peaks.by.genome.shuffled/${genome_id}.bed
  sort -k1,1 -k2,2n -o output/peaks.by.genome.shuffled/${genome_id}.bed output/peaks.by.genome.shuffled/${genome_id}.bed 
  scripts/annotate-intervals.py --flank 200 --offset 100 -g genomes/bed/${genome_id}.bed -b output/peaks.by.genome.shuffled/${genome_id}.bed -o output/peaks.by.genome.shuffled.annotated/${genome_id}.bed -c genomes/fasta/${genome_id}.fa.fai
done
