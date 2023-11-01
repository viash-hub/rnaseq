#!/bin/bash

sample="SRR6357070"
CURR=`pwd` 

echo ">>> Testing RSeQC bamstat module functionality"

viash ns build --setup cb --parallel -q bamstat

echo ">> Docker"
target/docker/rseqc/rseqc_bamstat/rseqc_bamstat \
    --input "$CURR/testData/test/$sample.bam" \
    --output output/bamstat

echo "> Checking whether output dir exists"
[[ ! -d output/bamstat ]] && echo "Output dir could not be found!" && exit 1

echo "> Checking if correct file is present"
[[ ! -f output/bamstat/"$sample".bamstat.txt ]] && echo "Bamstat summary file missing" && exit 1

echo ">> Nextflow"
nextflow run target/nextflow/rseqc/rseqc_bamstat/main.nf \
    --input "$CURR/testData/test/$sample.bam" \
    --output bamstat \
    --publish_dir nxf_output

echo "> Checking whether output dir exists"
[[ ! -d nxf_output/bamstat ]] && echo "Output dir could not be found!" && exit 1

echo "> Checking if correct file is present"
[[ ! -f nxf_output/bamstat/"$sample".bamstat.txt ]] && echo "Bamstat summary file missing" && exit 1
