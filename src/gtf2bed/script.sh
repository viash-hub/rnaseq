#!/bin/bash

set -eo pipefail 

filename=$(basename -- $par_gtf/*.gtf)
mkdir -p $par_bed_output

"$meta_resources_dir/gtf2bed" $par_gtf/$filename > $par_bed_output/${filename%%.*}.bed 