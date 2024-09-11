#!/bin/bash

set -eo pipefail

function clean_up {
  rm -rf "$tmpdir"
}
trap clean_up EXIT

tmpdir=$(mktemp -d "$meta_temp_dir/$meta_name-XXXXXXXX")

IFS="," read -ra input <<< $par_input
count=${#input[@]}

if $par_paired; then
  echo "Paired - $count"
  if [ $count -ne 2 ]; then
    echo "Paired end input requires two files"
    exit 1
  fi
else
  echo "Not Paired - $count"
  if [ $count -ne 1 ]; then
    echo "Single end input requires one file"
    exit 1
  fi
fi

fastqc -o $tmpdir ${input[*]} 

file1=$(basename -- "${input[0]}")
read1="${file1%.fastq*}"
file2=$(basename -- "${input[1]}")
read2="${file2%.fastq*}"

[[ -e "${tmpdir}/${read1}_fastqc.html" ]] && cp "${tmpdir}/${read1}_fastqc.html" $par_fastqc_html_1
[[ -e "${tmpdir}/${read2}_fastqc.html" ]] && cp "${tmpdir}/${read2}_fastqc.html" $par_fastqc_html_2
[[ -e "${tmpdir}/${read1}_fastqc.zip" ]] && cp "${tmpdir}/${read1}_fastqc.zip" $par_fastqc_zip_1
[[ -e "${tmpdir}/${read2}_fastqc.zip" ]] && cp "${tmpdir}/${read2}_fastqc.zip" $par_fastqc_zip_2
