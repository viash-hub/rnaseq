#!/bin/bash

echo ">>> Testing $meta_functionality_name"

echo ">>> Testing for paired-end reads"
"$meta_executable" \
    --paired true \
    --input $meta_resources_dir/SRR6357070_1.fastq.gz,$meta_resources_dir/SRR6357070_2.fastq.gz \
    --trim_html_1 SRR6357070_1.trimmed.html \
    --trim_html_2 SRR6357070_2.trimmed.html \
    --trim_zip_1 SRR6357070_1.trimmed.zip \
    --trim_zip_2 SRR6357070_2.trimmed.zip \
    --fastq_1 SRR6357070_1.trimmed.fastq.gz \
    --fastq_2 SRR6357070_2.trimmed.fastq.gz \
    --trim_log_1 SRR6357070_1.trimming_report.txt \
    --trim_log_2 SRR6357070_2.trimming_report.txt

echo ">> Checking if the correct files are present"
[[ ! -f SRR6357070_1.trimmed.html ]] || [[ ! -f SRR6357070_2.trimmed.html ]] && echo "Report file missing" && exit 1
[[ ! -f SRR6357070_1.trimmed.zip ]] || [[ ! -f SRR6357070_2.trimmed.zip ]] && echo "Zip file missing" && exit 1
[[ ! -f SRR6357070_1.trimmed.fastq.gz ]] || [[ ! -f SRR6357070_2.trimmed.fastq.gz ]] && echo "Trimmed reads file missing" && exit 1
[[ ! -f SRR6357070_2.trimming_report.txt ]] || [[ ! -f SRR6357070_2.trimming_report.txt ]] && echo "Trimming report log file missing" && exit 1


echo ">>> Testing for single-end reads"
"$meta_executable" \
    --paired false \
    --input $meta_resources_dir/SRR6357070_1.fastq.gz \
    --trim_html_1 SRR6357070_1.trimmed.html \
    --trim_zip_1 SRR6357070_1.trimmed.zip \
    --fastq_1 SRR6357070_1.trimmed.fastq.gz \
    --trim_log_1 SRR6357070_1.trimming_report.txt 

echo ">> Checking if the correct files are present"
[ ! -f SRR6357070_1.trimmed.html ] && echo "Report file missing" && exit 1
[ ! -f SRR6357070_1.trimmed.zip ] && echo "Zip file missing" && exit 1
[ ! -f SRR6357070_1.trimmed.fastq.gz ] && echo "Trimmed reads file missing" && exit 1
[ ! -f SRR6357070_2.trimming_report.txt ] && echo "Trimming report log file missing" && exit 1

echo ">>> Test finished successfully"
exit 0
