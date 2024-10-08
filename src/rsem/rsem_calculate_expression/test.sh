#!/bin/bash

echo ">>> Testing $meta_functionality_name"

tar -xavf $meta_resources_dir/rsem.tar.gz

echo ">>> Calculating expression"
"$meta_executable" \
  --id WT_REP1 \
  --strandedness reverse \
  --paired true \
  --input "$meta_resources_dir/SRR6357070_1.fastq.gz,$meta_resources_dir/SRR6357070_2.fastq.gz" \
  --index rsem \
  --extra_args "--star --star-output-genome-bam --star-gzipped-read-file --estimate-rspd --seed 1" \
  --counts_gene WT_REP1.genes.results \
  --counts_transctips WT_REP1.isoforms.results \
  --logs WT_REP1.log 
  
echo ">>> Checking whether output exists"
[ ! -f "WT_REP1.genes.results" ] && echo "Gene level expression counts file does not exist!" && exit 1
[ ! -s "WT_REP1.genes.results" ] && echo "Gene level expression counts file is empty!" && exit 1
[ ! -f "WT_REP1.log" ] && echo "Log file does not exist!" && exit 1
[ ! -s "WT_REP1.log" ] && echo "Log file is empty!" && exit 1

echo "All tests succeeded!"
exit 0
