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
html1=$(find $tmpdir/ -iname "*.read_1*.html")
html2=$(find $tmpdir/ -iname "*.read2*.html")
zip1=$(find $tmpdir/ -iname "*.read_1*.zip")
zip2=$(find $tmpdir/ -iname "*.read_2*.zip")

if [ -e "$html1" ]; then 
  cp $html1 $par_fastqc_html_1
fi

if [ -e "$html2" ]; then 
  cp $html2 $par_fastqc_html_2
fi

if [ -e "$zip1" ]; then
  cp $zip1 $par_fastqc_zip_1
fi

if [ -e "$zip2" ]; then
  cp $zip2 $par_fastqc_zip_2
fi

# # Version
# text="${meta_functionality_name}:
#   fastqc: $( fastqc --version | sed -e 's/FastQC v//g' )"

# if [ -e "$par_versions" ]; then
#     echo "$text" >> "$par_versions"
#     mv "$par_versions" "$par_updated_versions"
# else
#     echo "$text" > "$par_updated_versions"
# fi