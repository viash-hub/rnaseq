#!/bin/bash

set -eo pipefail

function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

tmpdir=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXXXX")

if [ $par_strandedness == 'forward' ]; then
    strandedness='--strandedness forward'
elif [ $par_strandedness == 'reverse' ]; then
    strandedness='--strandedness reverse'
else
    strandedness=''
fi

IFS="," read -ra input <<< $par_input

INDEX=`find -L $meta_resources_dir/ -name "*.grp" | sed 's/\.grp$//'`

rsem-calculate-expression \
    ${meta_cpus:+--num-theads $meta_cpus} \
    $strandedness \
    ${par_paired:+--paired-end} \
    $par_extra_args \
    ${input[*]} \
    $INDEX \
    $par_id
    
