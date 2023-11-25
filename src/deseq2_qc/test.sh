#!/bin/bash

# define input and output for script

counts="$meta_resources_dir/salmon.merged.gene_counts_length_scaled.tsv"

deseq2_output="deseq2"
pca_multiqc="star_salmon.pca.vals_mqc.tsv"
dists_multiqc="star_salmon.sample.dists_mqc.tsv"
extra_args="--id_col 1 --sample_suffix '' -outprefix deseq2 --count_col 3"
extra_args2="star_salmon"

# create temporary directory
tmpdir=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

# Run executable
echo "> Running $meta_functionality_name, writing to tmpdir $tmpdir."

"$meta_executable" \
    --count_file $counts \
    $extra_args \
    --cores 1 \
    --outdir $deseq2_output

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> asserting output has been created for paired read input"

[[ ! -d "$tmpdir/$deseq2_output" ]] && echo "$deseq2_output was not created" && exit 1
[[ ! -f "$tmpdir/$pca_multiqc" ]] && echo "$pca_multiqc was not created" && exit 1
[[ ! -f "$tmpdir/$dists_multiqc" ]] && echo "$dists_multiqc was not created" && exit 1

exit 0