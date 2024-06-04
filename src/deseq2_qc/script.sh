#!/bin/sh

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

Rscript "$meta_resources_dir/deseq2_qc.r" \
    --count_file $par_counts \
    --outdir $par_deseq2_output \
    --cores ${meta_cpus:-1} \
    $par_extra_args

if [ -f "$par_deseq2_output/R_sessionInfo.log" ]; then
    sed "s/deseq2_pca/${label_lower}_deseq2_pca/g" < $par_pca_header_multiqc > tmp.txt
    sed -i -e "s/DESeq2 PCA/${label_upper} DESeq2 PCA/g" tmp.txt
    cat tmp.txt $par_deseq2_output/*.pca.vals.txt > $par_pca_multiqc

    sed "s/deseq2_clustering/${label_lower}_deseq2_clustering/g" < $par_clustering_header_multiqc > tmp.txt
    sed -i -e "s/DESeq2 sample/${label_upper} DESeq2 sample/g" tmp.txt
    cat tmp.txt $par_deseq2_output/*.sample.dists.txt > $par_dists_multiqc
fi

# Version
# deseq2_ver=$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
# text="${meta_functionality_name}:
#     r-base: $(echo $(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
#     bioconductor-deseq2: ${deseq2_ver}"
# if [ -e "$par_versions" ]; then
#     echo "$text" >> "$par_versions"
#     mv "$par_versions" "$par_updated_versions"
# else
#     echo "$text" > "$par_updated_versions"
# fi