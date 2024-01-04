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
read -r -d '' text <<- END_VERSIONS
"${meta_functionality_name}":
  python: \$(python --version | sed 's/Python //g')
END_VERSIONS

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
else
    echo "$text" > "$par_versions"
fi