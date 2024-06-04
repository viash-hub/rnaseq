#!/bin/bash

set -eo pipefail 

bam_file="$(basename -- $par_bam_input)"
bai_file="$(basename -- $par_bai_input)"

echo "$bam_file"
echo "$bai_file"

tmpdir=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}

cp $par_bam_input $tmpdir/$bam_file
cp $par_bai_input $tmpdir/$bai_file

tin.py \
    -i $tmpdir/$bam_file \
    -r $par_refgene \
    -c $par_minimum_coverage \
    -n $par_sample_size \
    -s $par_subtract_background

mv ${bam_file%.*}.summary.txt $par_output_tin_summary
mv ${bam_file%.*}.tin.xls $par_output_tin

clean_up

# Version
# text="${meta_functionality_name}:
#     rseqc: $(tin.py --version | sed -e 's/tin.py //g')"

# if [ -e "$par_versions" ]; then
#     echo "$text" >> "$par_versions"
#     mv "$par_versions" "$par_updated_versions"
# else
#     echo "$text" > "$par_updated_versions"
# fi