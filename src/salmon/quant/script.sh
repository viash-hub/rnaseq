#!/bin/bash

set -eo pipefail
par_paired=true
input="testData/paired_end_test/SRR6357070.transcriptome_deduped_index.output"
par_alignment_mode=true
par_transcript_fasta="testData/test_output/ref.prepare_genome.transcript_fasta_uncompressed.fasta"
par_gtf="testData/test_output/ref.prepare_genome.gtf_uncompressed.gtf"
par_extra_salmon_quant_args=""
par_output="salmon_quant_test"
meta_cpus=2

# IFS="," read -ra input <<< $par_input
# if $par_paired; then
#     input_reads="-1 ${input[0]} -2 ${input[1]}"
# else 
#     input_reads="-r ${input[0]}"
# fi

reference="--index $par_star_index"
if $par_alignment_mode; then
    reference="-t $par_transcript_fasta/*"
    input_reads="-a $input/*.bam"
fi

strandedness_opts=('A' 'U' 'SF' 'SR' 'IS' 'IU' 'ISF' 'ISR' 'OS' 'OU' 'OSF' 'OSR' 'MS' 'MU' 'MSF' 'MSR')
strandedness='A'
# if $par_lib_type; then
#     if [ ${strandedness_opts[@]} =~ $par_lib_type ]; then
#         strandedness='$par_lib_type'
#     else 
#         echo "[Salmon Quant] Invalid library type specified '--libType=$lib_type', defaulting to auto-detection with '--libType=A'."
#     fi
# else
#     if $par_paired; then
#         strandedness='IU'
#     else 
#         strandedness='U'
#     fi
#     if [ $par_strandedness == 'forward' ]; then
#         if $par_paired; then
#             strandedness='ISF'
#         else 
#             strandedness='SF'
#         fi
#     else if [ $par_strandedness == 'reverse' ]; then 
#         if $par_paired; then
#             strandedness='ISR'
#         else 
#             strandedness='SR'
#         fi
#     fi
# fi

salmon quant --geneMap $par_gtf/* --threads $meta_cpus --libType=$strandedness $reference $input_reads $par_extra_salmon_quant_args -o $par_output

# if [ -f $prefix/aux_info/meta_info.json ]; then
#     cp $prefix/aux_info/meta_info.json "${prefix}_meta_info.json"
# fi