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
        cp $read1 $par_fastq_1
        cp $read2 $par_fastq_2
    fi
else
    echo "Not Paired - $read_count"
    if [ "$read_count" -ne 1 ]; then
        echo "Single end input requires one read file and one UMI pattern"
        exit 1
    else
        read1="$(basename -- ${input[0]})"
        trim_galore $par_extra_trimgalore_args --gzip --fastqc ${input[0]} -o $tmpdir
        read1=$(find $tmpdir/ -iname *trimmed.fq.gz*)
        cp $read1 $par_fastq_1
    fi
fi

mkdir -p $par_trim_log
mkdir -p $par_trim_html
mkdir -p $par_trim_zip
log=$(find $tmpdir/ -iname *report.txt)
html=$(find $tmpdir/ -iname *.html)
zip=$(find $tmpdir/ -iname *.zip)
cp -r $log $par_trim_log
cp -r $html $par_trim_html
cp -r $zip $par_trim_zip

# get total number of reads from the log ooutput file after trimming
# if [ $(ls -l $par_trim_log | wc -l) > 0 ]; then
#     if $par_paired; then
#         log_file=`grep -rl "shorter than the length cutoff" $par_trim_log`
#     else
#         log_file=$par_trim_log/*
#     fi
#     total_reads=`grep "[^0-9]* sequences processed in total" $log_file | sed 's/[^0-9]*//g'`
#     filtered_reads=`grep "shorter than the length cutoff" $log_file | sed -n 's/.*: \([0-9]\+\) .*/\1/p'`
#     par_trim_read_count=$(($total_reads - $filtered_reads)) 
# else
#     par_trim_read_count=$(($par_min_trimmed_reads + 1))
# fi
