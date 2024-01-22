echo "> Running $meta_functionality_name."

# define input and output for script
input_bam="$meta_resources_dir/qualimap_test.bam"
input_gtf="$meta_resources_dir/qualimap_test_annot.gtf"
output_dir="qualimap_output"
mkdir -p $output_dir

"$meta_executable" \
    --input "$input_bam" \
    --gtf "$input_gtf" \
    --output_dir "$output_dir"

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> Checking whether output dir and files exists"

[ ! -d "$output_dir" ] && echo "Output dir could not be found!" && exit 1
[ ! -d "$output_dir/raw_data_qualimapReport" ] && echo "Raw data folder could not be found!" && exit 1
[ -z $(ls -A "$output_dir/raw_data_qualimapReport") ] && echo "Raw data folder is missing output files" && exit 1
[ ! -f "$output_dir/qualimapReport.html" ] && echo "Qualimap report was not found" && exit 1
[ ! -s "$output_dir/qualimapReport.html" ] && echo "Qualimap report is empty" && exit 1

exit 0