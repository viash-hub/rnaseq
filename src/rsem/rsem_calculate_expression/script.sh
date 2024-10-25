#!/bin/bash

set -eo pipefail

function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

tmpdir=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXXXX")

[[ "$par_paired" == "false" ]] && unset par_paired

if [ $par_strandedness == 'forward' ]; then
    strandedness='--strandedness forward'
elif [ $par_strandedness == 'reverse' ]; then
    strandedness='--strandedness reverse'
else
    strandedness=''
fi

IFS=";" read -ra input <<< $par_input

INDEX=`find -L $par_index/ -name "*.grp" | sed 's/\.grp$//'`

rsem-calculate-expression \
    ${meta_cpus:+--num-threads $meta_cpus} \
    $strandedness \
    ${par_paired:+--paired-end} \
    $par_extra_args \
    ${input[*]} \
    $INDEX \
    $par_id
    
[[ -e "${par_id}.genes.results" ]] && mv "${par_id}.genes.results" $par_counts_gene
[[ -e "${par_id}id.isoforms.results" ]] && mv "${par_id}id.isoforms.results" $par_counts_transcripts
[[ -e "${par_id}.stat" ]] && mv "${par_id}.stat" $par_stat
# [[ -e "${par_id}.log" ]] && mv "${par_id}.log" $par_logs
[[ -e "${par_id}.STAR.genome.bam" ]] && mv "${par_id}.STAR.genome.bam" $par_bam_star
[[ -e "${par_id}.genome.bam" ]] && mv "${par_id}.genome.bam" $par_bam_genome
[[ -e "${par_id}.transcript.bam" ]] && mv "${par_id}.transcript.bam" $par_bam_transcript