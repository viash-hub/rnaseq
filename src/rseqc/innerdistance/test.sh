#!/bin/bash


## VIASH START
meta_executable="target/docker/rseqc/rseqc_inferexperiment/rseqc_inferexperiment"
meta_resources_dir="testData/test"
meta_temp_dir="output/inferexperiment"
## VIASH END

echo "> Running $meta_functionality_name for unpaired reads."

sample="SRR6357070"
genome="genome_gfp"

"$meta_executable" \
    --input "$meta_resources_dir/test/$sample.bam" \
    --refgene "$meta_resources_dir/test/$genome.bed" \
    --paired "false"

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> asserting no output has been created for unpaired read input"

[[ -f "$meta_temp_dir/$sample.inner_distance_stats.txt" ]] && echo "inner_distance_stats.txt has been created" && exit 1
[[ -f "$meta_temp_dir/$sample.inner_distance.txt" ]] && echo "inner_distance.txt has been created" && exit 1
[[ -f "$meta_temp_dir/$sample.inner_distance_freq.txt" ]] && echo "inner_distance_freq.txt has been created" && exit 1
[[ -f "$meta_temp_dir/$sample.inner_distance_plot.pdf" ]] && echo "inner_distance_plot.pdf has been created" && exit 1
[[ -f "$meta_temp_dir/$sample.inner_distance_plot.r" ]] && echo "inner_distance_plot.r has been created" && exit 1

echo "> Running $meta_functionality_name for unpaired reads."

sample="Pairend_nonStrandSpecific_36mer_Human_hg19"
genome="hg19_rRna"

"$meta_executable" \
    --input "$meta_resources_dir/test/$sample.bam" \
    --refgene "$meta_resources_dir/test/$genome.bed" \
    --paired "true"

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> asserting output has been created for paired read input"

[[! -f "$meta_temp_dir/$sample.inner_distance_stats.txt" ]] && echo "inner_distance_stats.txt was not created" && exit 1
[[! -f "$meta_temp_dir/$sample.inner_distance.txt" ]] && echo "inner_distance.txt was not reated" && exit 1
[[! -f "$meta_temp_dir/$sample.inner_distance_freq.txt" ]] && echo "inner_distance_freq.txt was not created" && exit 1
[[! -f "$meta_temp_dir/$sample.inner_distance_plot.pdf" ]] && echo "inner_distance_plot.pdf was not created" && exit 1
[[! -f "$meta_temp_dir/$sample.inner_distance_plot.r" ]] && echo "inner_distance_plot.r was not created" && exit 1

exit 0