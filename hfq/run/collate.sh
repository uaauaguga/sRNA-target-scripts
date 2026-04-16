#!/bin/bash
#metadata=data/PRJEB10931-ERP012236.metadata.txt
#dataset=GSE198671-SRP364136
#dataset=GSE216133-SRP403529
#dataset=GSE234792-PRJNA983241
#dataset=PRJEB10931-ERP012236
#dataset=GSE131520-SRP198991
dataset=GSE243246-PRJNA1017475
metadata=data/${dataset}.metadata.txt
mkdir output/bed.combined/$dataset
while read run_id label;do
  #echo $run_id $label
  [ -s output/bed/$dataset/${run_id}.bed ] || echo $run_id
  cat output/bed/$dataset/${run_id}.bed >> output/bed.combined/$dataset/${label}.bed
done < $metadata
