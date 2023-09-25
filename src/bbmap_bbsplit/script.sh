# #!/bin/bash

set -eo pipefail

avail_mem=3072

other_refs=()
while IFS="," read -r name path 
do
    other_refs+=("ref_$name=$path")
done < $par_bbsplit_fasta_list


if $par_only_build_index; then
    if [ -f "$par_primary_ref" ] && [ ${#other_refs[@]} -gt 0 ]; then
        bbsplit.sh -Xmx${avail_mem}M ref_primary=$par_primary_ref ${other_refs[@]} path=$par_bbsplit_index threads=$meta_cpus
    else
        echo "ERROR: Please specify as input a primary fasta file along with names and paths to non-primary fasta files."
    fi
else
    # IFS="," read -ra input <<< "$par_input"
    index_files=''
    if [ -e "$par_bbsplit_index" ]; then
        index_files="path=$par_bbsplit_index"
    elif [ -f "$par_primary_ref" ] && [ ${#other_refs[@]} -gt 0 ]; then
        index_files="ref_primary=$primary_ref ${other_refs[@]}"
    else
        echo "ERROR: Please either specify a BBSplit index as input or a primary fasta file along with names and paths to non-primary fasta files."
    fi
    # fastq_in=`$par_paired ? "in=${input[0]} in2=${input[1]}" : "in=${input}"`
    # fastq_out=`$par_paired ? "basename=${par_id}_%_#.fastq.gz" : "basename=${par_id}_%.fastq.gz"`
    read1=$(find $par_input -iname *_1_trimmed.fq.gz*)
    read2=$(find $par_input -iname *_2_trimmed.fq.gz*)
    if $par_paired; then
        fastq_in="in=$read1 in2=$read2"
        fastq_out="basename=${par_id}_%_#.fastq.gz"
    else
        fastq_in="in=$read1"
        fastq_out="basename=${par_id}_%.fastq.gz"
    fi

    bbsplit.sh -Xmx${avail_mem}M $index_files threads=$meta_cpus $fastq_in $fastq_out refstats=${par_id}.stats.txt
fi
cp "${par_id}_primary.fastq.gz" $par_filtered_output