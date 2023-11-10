#!/bin/bash

# define input and output for script
input_bam_1_unpaired="SRR6357070.bam"
input_bed_1="genome_gfp.bed"

input_bam_2_paired="Pairend_nonStrandSpecific_36mer_Human_hg19.bam"
input_bed_2="hg19_rRna.bed"

output_stats="inner_distance_stats.txt"
output_dist="inner_distance.txt"
output_plot="inner_distance_plot.pdf"
output_plot_r="inner_distance_plot.r"
output_freq="inner_distance_freq.txt"

# create temporary directory
tmpdir=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

# Run executable

echo "> Running $meta_functionality_name for unpaired reads, writing to tmpdir $tmpdir."

"$meta_executable" \
    --input "$meta_resources_dir/$input_bam_1_upaired" \
    --refgene "$meta_resources_dir/$input_bed_1" \
    --output_stats $tmpdir/$output_stats \
    --output_dist $tmpdir/$output_dist \
    --output_plot $tmpdir/$output_plot \
    --output_plot_r $tmpdir/$output_plot_r \
    --output_freq $tmpdir/$output_freq

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> asserting no output has been created for unpaired read input"

[[ -f $tmpdir/"$output_stats" ]] && echo "$output_stats has been created" && exit 1
[[ -f $tmpdir/"$output_dist" ]] && echo "$output_dist has been created" && exit 1
[[ -f $tmpdir/"$$output_plot" ]] && echo "$output_plot has been created" && exit 1
[[ -f $tmpdir/"$$output_plot_r" ]] && echo "$output_plot_r has been created" && exit 1
[[ -f $tmpdir/"$$soutput_freq" ]] && echo "$output_freq has been created" && exit 1

echo "> Running $meta_functionality_name for paired reads, writing to tmpdir $tmpdir."

"$meta_executable" \
    --input "$meta_resources_dir/$input_bam_2_paired" \
    --refgene "$meta_resources_dir/$input_bed_2" \
    --paired \
    --output_stats $tmpdir/$output_stats \
    --output_dist $tmpdir/$output_dist \
    --output_plot $tmpdir/$output_plot \
    --output_plot_r $tmpdir/$output_plot_r \
    --output_freq $tmpdir/$output_freq

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> asserting output has been created for paired read input"

[[ ! -f "$tmpdir/$output_stats" ]] && echo "$output_stats was not created" && exit 1
[[ ! -f "$tmpdir/$output_dist" ]] && echo "$output_dist was not created" && exit 1
[[ ! -f "$tmpdir/$output_plot" ]] && echo "$output_plot was not created" && exit 1
[[ ! -f "$tmpdir/$output_plot_r" ]] && echo "$output_plot_r was not created" && exit 1
[[ ! -f "$tmpdir/$output_freq" ]] && echo "$output_freq was not created" && exit 1

exit 0