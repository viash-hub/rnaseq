#!/bin/bash

echo "Testing $meta_functionality_name"
"$meta_executable" \
  --fasta "$meta_resources_dir/genome.fasta" \
  --sizes genome.fasta.sizes \
  --fai genome.fasta.fai

echo ">>> Checking whether output exists"
[ ! -f "genome.fasta.sizes" ] && echo "Chromosome lengths file does not exist!" && exit 1
[ ! -s "genome.fasta.sizes" ] && echo "Chromosome lengths file is empty!" && exit 1
[ ! -f "genome.fasta.fai" ] && echo "FASTA index file does not exist!" && exit 1
[ ! -s "genome.fasta.fai" ] && echo "FASTA index file does is empty!" && exit 1

echo "All tests succeeded!"
exit 0