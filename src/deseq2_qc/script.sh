#!/bin/sh

# # VIASH START
# par_counts="work/4e/709a2cba757cd071f1f7a7995d83c4/_viash_par/counts_gene_length_scaled_1/salmon.merged.gene_counts_length_scaled.tsv"
# par_pca_header_multiqc="assets/multiqc/deseq2_pca_header.txt"
# par_clustering_header_multiqc="assets/multiqc/deseq2_clustering_header.txt"
# meta_cpus=2
# par_extra_args2="star_salmon"
# par_extra_args="--id_col 1 --sample_suffix '' --outprefix deseq2 --count_col 3"
# par_deseq2_vst=true
# par_deseq2_output="deseq2_test"
# par_pca_multiqc="test_multiqc"
# par_dists_multiqc="dists_multiqc"
# # VIASH END

set -eo pipefail

if $par_deseq2_vst; then 
    par_extra_args+=" --vst TRUE"
fi

tolower() {
    case $1 in
        *[[:upper:]]*)
            printf "%s\n" "$1" | tr '[:upper:]' '[:lower:]'
            ;;
        *)
            printf "%s\n" "$1"
            ;;
    esac
}

toupper() {
    case $1 in
        *[[:lower:]]*)
            printf "%s\n" "$1" | tr '[:lower:]' '[:upper:]'
            ;;
        *)
            printf "%s\n" "$1"
            ;;
    esac
}

label_lower=$(tolower "$par_extra_args2")
label_upper=$(toupper "$par_extra_args2")

Rscript "$meta_resources_dir/deseq2_qc.r" --count_file $par_counts --outdir $par_deseq2_output --cores $meta_cpus $par_extra_args

if [ -f "$par_deseq2_output/R_sessionInfo.log" ]; then
    sed "s/deseq2_pca/${label_lower}_deseq2_pca/g" < $par_pca_header_multiqc > tmp.txt
    sed -i -e "s/DESeq2 PCA/${label_upper} DESeq2 PCA/g" tmp.txt
    cat tmp.txt $par_deseq2_output/*.pca.vals.txt > $par_pca_multiqc

    sed "s/deseq2_clustering/${label_lower}_deseq2_clustering/g" < $par_clustering_header_multiqc > tmp.txt
    sed -i -e "s/DESeq2 sample/${label_upper} DESeq2 sample/g" tmp.txt
    cat tmp.txt $par_deseq2_output/*.sample.dists.txt > $par_dists_multiqc
fi
