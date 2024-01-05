#!/bin/bash

set -eo pipefail

function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT 

tmpdir=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXXXX")

IFS="," read -ra input <<< "$par_input"
read_count="${#input[@]}"

if [ "$par_paired" == "true" ]; then
    echo "Paired - $read_count"
    if [ "$read_count" -ne 2 ]; then
        echo "Paired end input requires two read files and two UMI patterns"
        exit 1
    else
        trim_galore $par_extra_trimgalore_args --paired --gzip --fastqc ${input[0]} ${input[1]} -o $tmpdir 
        read1=$(find $tmpdir/ -iname *_1*.fq.gz*)
        read2=$(find $tmpdir/ -iname *_2*.fq.gz*)
        log1=$(find $tmpdir/ -iname *read1*report.txt)
        log2=$(find $tmpdir/ -iname *read2*report.txt)
        html1=$(find $tmpdir/ -iname *_1*.html)
        html2=$(find $tmpdir/ -iname *_2*.html)
        zip1=$(find $tmpdir/ -iname *_1*.zip)
        zip2=$(find $tmpdir/ -iname *_2*.zip)
        echo "$log1 - $log2"
        cp $read1 $par_fastq_1
        cp $read2 $par_fastq_2
        cp $log1 $par_trim_log_1
        cp $log2 $par_trim_log_2
        cp $html1 $par_trim_html_1
        cp $html2 $par_trim_html_2
        cp $zip1 $par_trim_zip_1
        cp $zip2 $par_trim_zip_2
    fi
else
    echo "Not Paired - $read_count"
    if [ "$read_count" -ne 1 ]; then
        echo "Single end input requires one read file and one UMI pattern"
        exit 1
    else
        read1="$(basename -- ${input[0]})"
        trim_galore $par_extra_trimgalore_args --gzip --fastqc ${input[0]} -o $tmpdir
        read=$(find $tmpdir/ -iname *trimmed.fq.gz*)
        log=$(find $tmpdir/ -iname *report.txt)
        html=$(find $tmpdir/ -iname *.html)
        zip=$(find $tmpdir/ -iname *.zip)
        cp $read $par_fastq_1
        cp $log $par_trim_log_1
        cp $html $par_trim_html_1
        cp $zip $par_trim_zip_1
    fi
fi

# Version
text="${meta_functionality_name}:
    trimgalore: $(echo $(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
    cutadapt: $(cutadapt --version)"

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
    mv "$par_versions" "$par_updated_versions"
else
    echo "$text" > "$par_updated_versions"
fi