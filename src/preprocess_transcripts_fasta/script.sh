#!/bin/bash

set -eo pipefail

command="cat"
if [ "${par_transcript_fasta##*.}"==".gz" ]; then
    command="zcat"
fi

$command $par_transcript_fasta | cut -d "|" -f1 > $par_fixed_fasta