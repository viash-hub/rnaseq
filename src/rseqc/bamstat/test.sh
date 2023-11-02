#!/bin/bash


## VIASH START
meta_executable="target/docker/rseqc/rseqc_bamstat/rseqc_bamstat"
meta_resources_dir="testData/test"
meta_temp_dir="output/bamstat"
## VIASH END

sample="SRR6357070"


"$meta_executable" \
    --input "$meta_resources_dir/$sample.bam" \
    --output "$meta_temp_dir"


[[ ! -d "$meta_temp_dir" ]] && echo "Output dir could not be found!" && exit 1

[[ ! -f "$meta_temp_dir/$sample".bamstat.txt ]] && echo "Bamstat summary file missing" && exit 1

exit 0