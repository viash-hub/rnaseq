echo "> Running $meta_functionality_name."

sample="SRR6357070"
genome="genome_gfp"

output_duplication_rate_plot_r="duplication_rate_plot.r"
output_duplication_rate_plot="duplication_rate_plot.pdf"

"$meta_executable" \
    --input "$meta_resources_dir/test/$sample.bam" \
    --output_duplication_rate_plot_r "$meta_temp_dir/$output_duplication_rate_plot_r" \
    --output_duplication_rate_plot "$meta_temp_dir/$output_duplication_rate_plot"

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> asserting  all output files were created"

[[ ! -f "$meta_temp_dir/$output_duplication_rate_plot_r" ]] && echo "$output_duplication_rate_plot_r was not created" && exit 1
[[ ! -f "$meta_temp_dir/$output_duplication_rate_plot" ]] && echo "$output_duplication_rate_plot was not created" && exit 1

exit 0