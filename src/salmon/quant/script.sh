#!/bin/bash

set -eo pipefail

if $par_alignment_mode; then
    reference="-t $par_transcript_fasta"
    input_reads="-a $par_input"
else
    IFS="," read -ra input <<< $par_input
    if $par_paired; then
        input_reads="-1 ${input[0]} -2 ${input[1]}"
    else 
        input_reads="-r ${input[0]}"
    fi
    reference="--index $par_salmon_index"
fi

strandedness_opts=('A' 'U' 'SF' 'SR' 'IS' 'IU' 'ISF' 'ISR' 'OS' 'OU' 'OSF' 'OSR' 'MS' 'MU' 'MSF' 'MSR')
strandedness='A'
if [ $par_lib_type != '' ]; then
    if [[ ${strandedness_opts[@]} =~ $par_lib_type ]]; then
        strandedness=$par_lib_type
    else 
        echo "[Salmon Quant] Invalid library type specified '--libType=$lib_type', defaulting to auto-detection with '--libType=A'."
    fi
else
    if [ $par_strandedness == 'forward' ]; then
        if $par_paired; then
            strandedness='ISF'
        else 
            strandedness='SF'
        fi
    elif [ $par_strandedness == 'reverse' ]; then 
        if $par_paired; then
            strandedness='ISR'
        else 
            strandedness='SR'
        fi
    else 
        if $par_paired; then
            strandedness='IU'
        else 
            strandedness='U'
        fi
    fi
fi

salmon quant \
    --geneMap $par_gtf \
    --threads $meta_cpus \
    --libType=$strandedness \
    $reference \
    $input_reads \
    $par_extra_salmon_quant_args \
    -o $par_output

if [ -f "$par_output/aux_info/meta_info.json" ]; then
    cp "$par_output/aux_info/meta_info.json" $par_json_info
fi

# Version
read -r -d '' text <<- END_VERSIONS
"${meta_functionality_name}":
    salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
END_VERSIONS

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
else
    echo "$text" > "$par_versions"
fi