#!/bin/bash

set -eo pipefail

function clean_up {
  rm -rf "$tmpdir"
}
trap clean_up EXIT

tmpdir=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXXXX")

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
html1=$(find $tmpdir/ -iname ${read1}_fastqc.html)
html2=$(find $tmpdir/ -iname ${read2}_fastqc.html)
zip1=$(find $tmpdir/ -iname ${read1}_fastqc.zip)
zip2=$(find $tmpdir/ -iname ${read2}_fastqc.zip)
[ -e "$html1" ] && cp $html1 $par_fastqc_html_1
[ -e "$html2" ] && cp $html2 $par_fastqc_html_2
[ -e "$zip1" ] && cp $zip1 $par_fastqc_zip_1
[ -e "$zip2" ] && cp $zip2 $par_fastqc_zip_2

