#!/bin/bash

# Run executable
echo "> Running $meta_functionality_name"

"$meta_executable" \
    --counts $meta_resources_dir/counts.tsv \
    --pca_header_multiqc $meta_resources_dir/deseq2_pca_header.txt \
    --clustering_header_multiqc $meta_resources_dir/deseq2_clustering_header.txt \
    --extra_args "--id_col 1 --sample_suffix '' --outprefix deseq2 --count_col 2" \
    --extra_args2 "test" \
    --deseq2_output "deseq2/" \
    --pca_multiqc pca.vals_mqc.tsv \
    --dists_multiqc sample.dists_mqc.tsv

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