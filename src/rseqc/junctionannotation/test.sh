#!/bin/bash

# define input and output for script
input_bam="SRR6357070.bam"
input_bed="genome_gfp.bed"

output_junction_bed="junction_annotation.bed"
output_junction_interact="junction_annotation.Interact.bed"
output_junction_sheet="junction_annotation.xls"
output_plot_r="unction_annotation_plot.r"
output_splice_events_plot="splice_events.pdf"
output_splice_junctions_plot="splice_junctions_plot.pdf"
output_log="junction_annotation.log"

# create temporary directory
tmpdir=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

# run executable and test
echo "> Running $meta_functionality_name, writing to tmpdir $tmpdir."

"$meta_executable" \
    --input "$meta_resources_dir/$input_bam" \
    --refgene "$meta_resources_dir/$input_bed" \
    --output_log "$tmpdir/$output_log" \
    --output_plot_r "$tmpdir/$output_plot_r" \
    --output_junction_bed "$tmpdir/$output_junction_bed" \
    --output_junction_interact "$tmpdir/$output_junction_interact" \
    --output_junction_sheet "$tmpdir/$output_junction_sheet" \
    --output_splice_events_plot "$tmpdir/$output_splice_events_plot" \
    --output_splice_junctions_plot "$tmpdir/$output_splice_junctions_plot" 

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> asserting  all output files were created"

[[ ! -f "$tmpdir/$output_log" ]] && echo "$output_log was not created" && exit 1
[[ ! -f "$tmpdir/$output_junction_bed" ]] && echo "$output_junction_bed was not created" && exit 1
[[ ! -f "$tmpdir/$output_junction_interact" ]] && echo "$output_junction_interact was not created" && exit 1
[[ ! -f "$tmpdir/$output_junction_sheet" ]] && echo "$output_junction_sheet was not created" && exit 1
[[ ! -f "$tmpdir/$output_plot_r" ]] && echo "$output_plot_r was not created" && exit 1
[[ ! -f "$tmpdir/$output_splice_events_plot" ]] && echo "$output_splice_events_plot was not created" && exit 1
[[ ! -f "$tmpdir/$output_splice_junctions_plot" ]] && echo "$output_splice_junctions_plot was not created" && exit 1

exit 0