#!/bin/bash
for genome_id in $(cat  genome-ids.txt);do
  ls output/peaks.by.genome/${genome_id}.sorted.merged.bed
  scripts/annotate-intervals.py --flank 200 --offset 100 -g genomes/bed/${genome_id}.bed -b output/peaks.by.genome/${genome_id}.sorted.merged.bed -o output/peaks.by.genome.annotated/${genome_id}.bed -c genomes/fasta/${genome_id}.fa.fai
done
