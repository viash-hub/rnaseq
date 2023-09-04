#!/bin/bash

set -eo pipefail

# par_input="/home/nirmayi/data_intuitive/rnaseq.vsh/testData/reference/hisat2.tar.gz"

prefix=`echo "${par_input%%.tar.gz}"`

mkdir -p $prefix

## Ensures --strip-components only applied when top level of tar contents is a directory
## If just files or multiple directories, place all in prefix
if [[ $(tar -taf ${par_input} | grep -o -P "^.*?\\/" | uniq | wc -l) -eq 1 ]]; then
    tar -C $prefix --strip-components 1 -xavf $par_input --no-same-owner
else
    tar -C $prefix -xavf $par_input --no-same-owner
fi