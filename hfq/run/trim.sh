#!/bin/bash
#dataset=GSE234792-PRJNA983241
dataset=GSE243246-PRJNA1017475
for run_id in $(cat data/${dataset}.txt);do
  bash scripts/trim.sh $dataset $run_id
done

