# #!/bin/bash

set -eo pipefail
# par_input="/home/nirmayi/data_intuitive/rnaseq.vsh/testData/single_end_test/SRR6357070_1.trimgalore.output/SRR6357070_1.umitools_extract.output_1_trimmed.fq.gz"
# par_bbsplit_fasta_list="/home/nirmayi/data_intuitive/rnaseq.vsh/testData/reference/bbsplit_fasta_list.txt"
# par_only_build_index=true
# par_primary_ref="/home/nirmayi/data_intuitive/rnaseq.vsh/testData/reference/genome.fasta"
# par_bbsplit_index="/home/nirmayi/data_intuitive/rnaseq.vsh/testData/bbsplit"
mkdir -p $par_bbsplit_index
meta_cpus=2

avail_mem=3072

other_refs=()
while IFS="," read -r name path 
do
    other_refs+=("ref_$name=$path")
done < $par_bbsplit_fasta_list

if [ $par_only_build_index ]; then
    # if [ (-f "$par_primary_ref") && (-f "$par_other_refs") ]; then
    echo "bbsplit.sh -Xmx${avail_mem}M ref_primary=$par_primary_ref ${other_refs[@]} path=$par_bbsplit_index threads=$meta_cpus"
    # else
        # echo "ERROR: Please specify as input a primary fasta file along with names and paths to non-primary fasta files."
    # fi
# else
#     IFS="," read -ra input <<< "$par_input"
#     index_files=''
#     if [ $par_bbsplit_index ]; then
#         index_files="path=$par_bbsplit_index"
#     elif [ $par_primary_ref && $par_other_refs ]; then
#         index_files="ref_primary=${$primary_ref} ${other_refs[@]}"
#     else
#         echo "ERROR: Please either specify a BBSplit index as input or a primary fasta file along with names and paths to non-primary fasta files."
#     fi
#     fastq_in=`$par_paired ? "in=${input[0]} in2=${input[1]}" : "in=${input}"`
#     fastq_out=`$par_paired ? "basename=${par_id}_%_#.fastq.gz" : "basename=${par_id}_%.fastq.gz"`

#     echo "bbsplit.sh -Xmx${avail_mem}M $index_files threads=$meta_cpus $fastq_in $fastq_out refstats=${par_id}.stats.txt"
fi