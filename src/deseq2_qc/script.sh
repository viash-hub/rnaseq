#!/bin/bash

# par_counts="testData/paired_end_test/SRR6357071.genome_alignment_and_quant.salmon_tximport/salmon.merged.gene_counts_length_scaled.tsv"
# par_pca_header_multiqc="assets/multiqc/deseq2_pca_header.txt"
# par_clustering_header_multiqc="assets/multiqc/deseq2_clustering_header.txt"
# meta_cpus=2
# par_extra_args2="star_salmon"
# par_extra_args="--id_col 1 --sample_suffix '' --outprefix deseq2 --count_col 3"
# par_deseq2_vst=true

set -eo pipefail

if $par_deseq2_vst; then 
    par_extra_args+=" --vst TRUE"
fi
label_lower=${par_extra_args2,,}
label_upper=${par_extra_args2^^}

"$meta_resources_dir/deseq2_qc.r" --count_file $par_counts --outdir $par_deseq2_output --cores $meta_cpus $par_extra_args

if [ -f "$par_deseq2_output/R_sessionInfo.log" ]; then
    sed "s/deseq2_pca/${label_lower}_deseq2_pca/g" < $par_pca_header_multiqc > tmp.txt
    sed -i -e "s/DESeq2 PCA/${label_upper} DESeq2 PCA/g" tmp.txt
    cat tmp.txt $par_deseq2_output/*.pca.vals.txt > $par_pca_multiqc

    sed "s/deseq2_clustering/${label_lower}_deseq2_clustering/g" < $par_clustering_header_multiqc > tmp.txt
    sed -i -e "s/DESeq2 sample/${label_upper} DESeq2 sample/g" tmp.txt
    cat tmp.txt $par_deseq2_output/*.sample.dists.txt > $par_dists_multiqc
fi
