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
html1=$(find $tmpdir/ -iname *read1*fastqc.html)
html2=$(find $tmpdir/ -iname *read2*fastqc.html)
zip1=$(find $tmpdir/ -iname *read1*fastqc.zip)
zip2=$(find $tmpdir/ -iname *read2*fastqc.zip)
cp $html1 $par_fastqc_html_1
cp $html2 $par_fastqc_html_2
cp $zip1 $par_fastqc_zip_1
cp $zip2 $par_fastqc_zip_2

