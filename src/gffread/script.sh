#!/bin/bash

set -eo pipefail

prefix=`echo "${par_input%%.gff*}"`
# par_output="$prefix.gtf"

gffread $par_input -o $par_output