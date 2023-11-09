echo "> Running $meta_functionality_name."

# define input and output for script

gunzip "$meta_resources_dir/genes.gtf.gz"

input_bam="$meta_resources_dir/SRR6357070.bam"
input_gtf="$meta_resources_dir/genes.gtf"

# create temporary directory
tmpdir=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

"$meta_executable" \
    --input "$input_bam" \
    --gtf "$input_gtf" \
    --output "$tmpdir/$output_dir"

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> Checking whether output dir and files exists"

[[ ! -d "$tmpdir/$output_dir" ]] && echo "Output dir could not be found!" && exit 1
[[ ! -d "$tmpdir/$output_dir/raw_data_qualimapReport" ]] && echo "Raw data folder could not be found!" && exit 1
[[ ! -f "$tmpdir/$output_dir/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt" ]] && echo "raw data folder is missing output files" && exit 1
[[ ! -f "$tmpdir/$output_dir/qualimapReport.html" ]] && echo "qualimap report was not found" && exit 1

exit 0