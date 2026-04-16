#!/bin/bash
#while read dataset asm_id;do
#  #echo $dataset $asm_id
#  cat output/peaks/${dataset}.bed >> output/peaks.by.genome/${asm_id}.bed
#done < dataset2genome.txt

for genome_id in $(cat dataset2genome.txt | cut -f 2 | sort -u);do
  echo $genome_id
  #[ -s output/peaks.by.genome/${genome_id}.sorted.bed ] || 
  cat output/peaks.by.genome/${genome_id}.bed | awk 'BEGIN{FS="\t";OFS="\t";}$7<0.0005{print}' | sort -k1,1 -k2,2n | cut -f 1-7 > output/peaks.by.genome/${genome_id}.sorted.bed
  bedtools merge -s -d 1 -c 4,5,6,7 -o count,max,distinct,min -i output/peaks.by.genome/${genome_id}.sorted.bed > output/peaks.by.genome/${genome_id}.sorted.merged.bed
done
