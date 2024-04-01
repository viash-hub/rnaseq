#!/bin/bash

echo ">>> Testing $meta_functionality_name"

echo ">>> Generating Kallisto index"
kallisto index \
    -i index \
    $meta_resources_dir/transcriptome.fasta

echo ">>> Testing for paired-end reads"
"$meta_executable" \
  --index index \
  --paired true \
  --strandedness reverse \
  --output paired_end_test \
  --input "SRR6357070_1.fastq.gz,SRR6357070_2.fastq.gz" \
  --log quant_pe.log \
  --run_info pe_run_info.json 

echo ">>> Checking whether output exists"
[ ! -d "paired_end_test" ] && echo "Kallisto results do not exist!" && exit 1
[ ! -f "quant_pe.log" ] && echo "quant_pe.log does not exist!" && exit 1
[ ! -s "quant_pe.log" ] && echo "quant_pe.log is empty!" && exit 1
[ ! -f "pe_run_info.json" ] && echo "pe_run_info.json does not exist!" && exit 1
[ ! -s "pe_run_info.json" ] && echo "pe_run_info.json is empty!" && exit 1
[ ! -f "paired_end_test/abundance.tsv" ] && echo "abundance.tsv does not exist!" && exit 1
[ ! -s "paired_end_test/abundance.tsv" ] && echo "abundance.tsv is empty!" && exit 1
[ ! -f "paired_end_test/abundance.h5" ] && echo "abundance.h5 does not exist!" && exit 1
[ ! -s "paired_end_test/abundance.h5" ] && echo "abundance.h5 is empty!" && exit 1

echo ">>> Testing for single-end reads"
"$meta_executable" \
  --index index \
  --paired false \
  --strandedness "reverse" \
  --output single_end_test \
  --input "SRR6357070_1.fastq.gz" \
  --log quant_se.log \
  --run_info se_run_info.json \
  --fragment_length 101 \
  --fragment_length_sd 50

echo ">>> Checking whether output exists"
[ ! -d "single_end_test" ] && echo "Kallisto results do not exist!" && exit 1
[ ! -f "quant_se.log" ] && echo "quant_se.log does not exist!" && exit 1
[ ! -s "quant_se.log" ] && echo "quant_se.log is empty!" && exit 1
[ ! -f "se_run_info.json" ] && echo "se_run_info.json does not exist!" && exit 1
[ ! -s "se_run_info.json" ] && echo "se_run_info.json is empty!" && exit 1
[ ! -f "single_end_test/abundance.tsv" ] && echo "abundance.tsv does not exist!" && exit 1
[ ! -s "single_end_test/abundance.tsv" ] && echo "abundance.tsv is empty!" && exit 1
[ ! -f "single_end_test/abundance.h5" ] && echo "abundance.h5 does not exist!" && exit 1
[ ! -s "single_end_test/abundance.h5" ] && echo "abundance.h5 is empty!" && exit 1

echo "All tests succeeded!"
exit 0
