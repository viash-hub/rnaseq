#!/bin/bash

set -eo pipefail

prefix=(echo "${par_input%%.tar.gz}")

mkdir -p $prefix

## Ensures --strip-components only applied when top level of tar contents is a directory
## If just files or multiple directories, place all in prefix
if [[ \$(tar -taf ${par_input} | grep -o -P "^.*?\\/" | uniq | wc -l) -eq 1 ]]; then
    tar -C $prefix --strip-components 1 -xavf $par_input \
        # args $args2
else
    tar -C $prefix -xavf $par_input \
        # $args $args2
fi