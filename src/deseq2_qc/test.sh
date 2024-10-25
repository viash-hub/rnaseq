#!/bin/bash

# Run executable
echo "> Running $meta_functionality_name"

"$meta_executable" \
    --counts $meta_resources_dir/counts.tsv \
    --id_col 1 \
    --sample_suffix '' \
    --outprefix deseq2 \
    --count_col 2 \
    --deseq2_output "deseq2/" \
    --pca_multiqc pca.vals_mqc.tsv \
    --sample_dists_multiqc sample.dists_mqc.tsv \
    --outdir deseq2

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> Check whether output exists"

[ ! -d "deseq2" ] && echo "deseq2 was not created" && exit 1
[ -z "$(ls -A 'deseq2')" ] && echo "deseq2 is empty" && exit 1
[ ! -f "pca.vals_mqc.tsv" ] && echo "pca.vals_mqc.tsv was not created" && exit 1
[ ! -s "pca.vals_mqc.tsv" ] && echo "pca.vals_mqc.tsv is empty" && exit 1
[ ! -f "sample.dists_mqc.tsv" ] && echo "sample.dists_mqc.tsv was not created" && exit 1
[ ! -s "sample.dists_mqc.tsv" ] && echo "sample.dists_mqc.tsv is empty" && exit 1

exit 0
