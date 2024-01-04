#!/bin/bash

set -eo pipefail

filename="$(basename -- "$par_transcript_fasta")"

if [ ${filename##*.} == "gz" ]; then
    zcat $par_transcript_fasta | cut -d "|" -f1 > $par_output
else 
    cat $par_transcript_fasta | cut -d "|" -f1 > $par_output
fi
