#!/bin/bash

## VIASH START
meta_resources_dir=/src/dupradar
meta_cpus=10
## VIASH END

set -exo pipefail 

filename="$(basename -- $par_input)"

function num_strandness {
    if [ $par_strandedness == 'none' ]; then echo 1
    elif [ $par_strandedness == 'forward' ]; then echo 1
    elif [ $par_strandedness == 'reverse' ]; then echo 2
    else echo "strandedness must be none, forward or reverse." && \
        exit 1
    fi
}

function str_paired {
    if $par_paired; then echo paired
    else echo single
    fi
}

Rscript $meta_resources_dir/dupRadar.r \
    $par_input \
    $filename \
    $par_gtf_annotation \
    $(num_strandness) \
    $par_paired \
    $par_cpus

mv "$filename"_dupMatrix.txt $par_output_dupmatrix
mv "$filename"_dup_intercept_mqc.txt $par_output_dup_intercept_mqc
mv "$filename"_duprateExpBoxplot.pdf $par_output_duprate_exp_boxplot
mv "$filename"_duprateExpDens.pdf $par_output_duprate_exp_densplot
mv "$filename"_duprateExpDensCurve_mqc.txt $par_output_duprate_exp_denscurve_mqc
mv "$filename"_expressionHist.pdf $par_output_expression_histogram
mv "$filename"_intercept_slope.txt $par_output_intercept_slope
