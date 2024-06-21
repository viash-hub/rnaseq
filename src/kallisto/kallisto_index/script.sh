#!/bin/bash

set -eo pipefail

kallisto index \
    ${par_pseudo_aligner_kmer_size:+-k $par_pseudo_aligner_kmer_size} \
    -i $par_kallisto_index \
    $par_transcriptome_fasta
