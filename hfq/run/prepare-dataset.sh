#!/bin/bash

for genome_id in $(cat training-genome-ids.txt);do
  cat output/binding.sites.by.genome/${genome_id}.fa >> dataset/training.fa
done

for genome_id in $(cat testing-genome-ids.txt);do
   cat output/binding.sites.by.genome/${genome_id}.fa >> dataset/testing.fa
done
