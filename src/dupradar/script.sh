#!/bin/bash

## VIASH START
meta_resources_dir=/src/dupradar
meta_cpus=10
## VIASH END

set -exo pipefail 

prefix=$(openssl rand -hex 8)

function num_strandness {
    if [ -z $par_strandedness ]; then echo 1
    elif [ $par_strandedness == 'forward' ]; then echo 1
    elif [ $par_strandedness == 'reverse' ]; then echo 2
    else echo "strandedness must be null, forward or reverse." && \
        exit 1
    fi
}

Rscript $meta_resources_dir/dupRadar.r \
    $par_input \
    $prefix \
    $par_gtf_annotation \
    $(num_strandness) \
    $par_paired \
    ${meta_cpus:-1}

mv "$prefix"_dupMatrix.txt $par_output_dupmatrix
mv "$prefix"_dup_intercept_mqc.txt $par_output_dup_intercept_mqc
mv "$prefix"_duprateExpBoxplot.pdf $par_output_duprate_exp_boxplot
mv "$prefix"_duprateExpDens.pdf $par_output_duprate_exp_densplot
mv "$prefix"_duprateExpDensCurve_mqc.txt $par_output_duprate_exp_denscurve_mqc
mv "$prefix"_expressionHist.pdf $par_output_expression_histogram
mv "$prefix"_intercept_slope.txt $par_output_intercept_slope
