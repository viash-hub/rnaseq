#!/bin/bash

set -eo pipefail

add_name="$(basename -- $par_additional_fasta)"

# Use fasta2gtf.py to generate a GTF annotation file from the FASTA file
python3 "$meta_resources_dir/fasta2gtf.py" \
  -o ${add_name%%.*}.gtf \
  ${par_biotype:+-b $par_biotype} \
  $biotype_name \
  $par_additional_fasta

cat $par_fasta $par_additional_fasta > $par_fasta_output
cat $par_gtf ${add_name%%.*}.gtf > $par_gtf_output

# Version

text="${meta_functionality_name}:
  python: $(python3 --version | sed 's/Python //g')"

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
    mv "$par_versions" "$par_updated_versions"
else
    echo "$text" > "$par_updated_versions"
fi