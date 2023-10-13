#!/bin/bash

set -eo pipefail

filename="$(basename -- $par_input/*.gff)"
mkdir -p $par_output

gffread $par_input -o $par_output/${filename%%.*}.gtf