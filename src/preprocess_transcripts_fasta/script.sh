#!/bin/bash

set -eo pipefail

filename="$(basename -- "$par_transcript_fasta")"
mkdir -p $par_output

if [ ${filename##*.} == "gz" ]; then
    zcat $par_transcript_fasta/$filename | cut -d "|" -f1 > $par_output/${filename%.*}
else 
    cat $par_transcript_fasta/$filename | cut -d "|" -f1 > $par_output/$filename
fi