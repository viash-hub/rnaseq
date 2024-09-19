#!/bin/bash

set -eo pipefail

IFS="," read -ra input <<< $par_input

single_end_params=''
if [ $par_paired == "false" ]; then
    if [[ $par_fragment_length < 0 ]] || [[ ! $fragment_length_sd < 0 ]]; then
        echo "fragment_length and fragment_length_sd must be set for single-end data"
        exit 1
    fi
    single_end_params="--single --fragment-length $par_fragment_length --sd $par_fragment_length_sd"
fi

strandedness=''
if [[ "$par_extra_args" != *"--fr-stranded"* ]] && [[ "$par_extra_args" != *"--rf-stranded"* ]]; then
    if [ "$par_strandedness" == 'forward' ]; then
        strandedness='--fr-stranded'
    elif [ "$par_strandedness" == 'reverse' ]; then
        strandedness='--rf-stranded'
    fi
fi

mkdir -p $par_output

kallisto quant \
    ${meta_cpus:+--threads $meta_cpus} \
    --index $par_index \
    ${par_gtf:+--gtf $par_gtf} \
    ${par_chromosomes:+--chromosomes $par_chromosomes} \
    $single_end_params \
    $strandedness \
    $par_extra_args \
    -o $par_output \
    ${input[*]} 2> >(tee -a ${par_output}/kallisto_quant.log >&2)

mv ${par_output}/kallisto_quant.log ${par_log}
mv ${par_output}/run_info.json ${par_run_info}
cp ${par_output}/abundance.tsv ${par_quant_results_file}
