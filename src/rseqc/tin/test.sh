echo "> Running $meta_functionality_name."

sample="SRR6357070"
genome="genome_gfp"

output_tin="tin.xls"
output_tin_summary="tin_summary.txt"

"$meta_executable" \
    --input "$meta_resources_dir/test/$sample.bam" \
    --refgene "$meta_resources_dir/test/$genome.bed" \
    --output_tin "$meta_temp_dir/$output_tin" \
    --output_tin_summary "$meta_temp_dir/$output_tin_summary"

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> asserting  all output files were created"

[[ ! -f "$meta_temp_dir/$output_tin" ]] && echo "$output_tin was not created" && exit 1
[[ ! -f "$meta_temp_dir/$output_tin_summary" ]] && echo "$output_tin_summary was not created" && exit 1

exit 0