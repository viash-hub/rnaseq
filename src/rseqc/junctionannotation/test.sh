#!/bin/bash

echo "> Running $meta_functionality_name."

sample="SRR6357070"
genome="genome_gfp"

"$meta_executable" \
    --input "$meta_resources_dir/test/$sample.bam" \
    --refgene "$meta_resources_dir/test/$genome.bed" \

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> asserting  all output files were created"

[[! -f "$meta_temp_dir/$output_log" ]] && echo "$output_log was not created" && exit 1
[[! -f "$meta_temp_dir/$par_output_junction_bed" ]] && echo "$par_output_junction_bed was not created" && exit 1
[[! -f "$meta_temp_dir/$par_output_junction_interact" ]] && echo "$par_output_junction_interact was not created" && exit 1
[[! -f "$meta_temp_dir/$par_output_junction_sheet" ]] && echo "$par_output_junction_sheet was not created" && exit 1
[[! -f "$meta_temp_dir/$par_output_plot_r" ]] && echo "$par_output_plot_r was not created" && exit 1
[[! -f "$meta_temp_dir/$par_output_splice_events_plot" ]] && echo "$par_output_splice_events_plot was not created" && exit 1
[[! -f "$meta_temp_dir/$par_output_splice_junctions_plot" ]] && echo "$par_output_splice_junctions_plot was not created" && exit 1

exit 0