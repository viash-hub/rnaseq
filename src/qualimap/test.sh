echo "> Running $meta_functionality_name."

sample="SRR6357070"
genome="genes"
output_dir="qualimap_output"

gunzip "$meta_resources_dir/reference/genes.gtf.gz"

"$meta_executable" \
    --input "$meta_resources_dir/test/$sample.bam" \
    --gtf "$meta_resources_dir/reference/$genome.gtf" \
    --output "$output_dir"

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> Checking whether output dir and files exists"

[[ ! -d "$output_dir" ]] && echo "Output dir could not be found!" && exit 1
[[ ! -d "$output_dir/raw_data_qualimapReport" ]] && echo "Raw data folder could not be found!" && exit 1
[[ ! -f "$output_dir/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt" ]] && echo "raw data folder is missing output files" && exit 1
[[ ! -f "$output_dir/qualimapReport.html" ]] && echo "qualimap report was not found" && exit 1

exit 0