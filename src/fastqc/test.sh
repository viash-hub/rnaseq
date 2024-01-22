#!/bin/bash

echo ">>> Testing $meta_functionality_name"

echo ">>> Testing for paired-end reads"

"$meta_executable" \
    --paired true \
    --input $meta_resources_dir/SRR6357070_1.fastq.gz,$meta_resources_dir/SRR6357070_2.fastq.gz \
    --fastqc_html_1  SRR6357070_1.html \
    --fastqc_html_2  SRR6357070_2.html \
    --fastqc_zip_1  SRR6357070_1.zip \
    --fastqc_zip_2  SRR6357070_2.zip

echo ">> Checking if the correct files are present"
[[ ! -f "SRR6357070_1.html" ]] || [[ ! -f "SRR6357070_2.html" ]] && echo "Report file missing" && exit 1
[[ ! -s "SRR6357070_1.html" ]] || [[ ! -s "SRR6357070_2.html" ]] && echo "Report file empty" && exit 1
[[ ! -f "SRR6357070_1.zip" ]] || [[ ! -f "SRR6357070_2.zip" ]] && echo "Zip file missing" && exit 1

rm SRR6357070_1.html SRR6357070_2.html SRR6357070_1.zip SRR6357070_2.zip

echo ">>> Testing for single-end reads"
"$meta_executable" \
    --paired false \
    --input $meta_resources_dir/SRR6357070_1.fastq.gz \
    --fastqc_html_1 SRR6357070_1.html \
    --fastqc_zip_1 SRR6357070_1.zip 
    
echo ">> Checking if the correct files are present"
[ ! -f "SRR6357070_1.html" ] && echo "Report file missing" && exit 1
[ ! -s "SRR6357070_1.html" ] && echo "Report file empty" && exit 1
[ ! -f "SRR6357070_1.zip" ] && echo "Zip file missing" && exit 1

echo ">>> Test finished successfully"
exit 0
