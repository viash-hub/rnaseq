#!/bin/bash

set -eo pipefail

par_input="testData/test/SRR6357070_1.fastq.gz"
par_output="umi_extract_test/"

mkdir -p $par_output

# if [ "$par_mode" == "dir" ]; then
#   par_input="$par_input/*.fastq.gz"
#   filename=`"echo "$par_input" | awk '{gsub(/.*[/]|[.].*/, "", $0)} 1'`
# fi

# single end
filename="${par_input%%.*}"
eval umi_tools extract -I "$par_input" -S "${filename}.umi_extract.fastq.gz" 

# paired end
