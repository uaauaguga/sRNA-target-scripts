#!/bin/bash
mkdir -p hfq-scores-0412
for genome_id in $(cat genome-ids.txt );do
  scripts/leader-scoring.py -f assemblies/${genome_id}.fa -b CDS/${genome_id}.bed -o hfq-scores-0412/${genome_id}.txt > hfq-scores-0412/${genome_id}.log 2>&1 
done
