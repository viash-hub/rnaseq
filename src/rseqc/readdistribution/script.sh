#!/bin/bash

set -eo pipefail 

read_distribution.py \
    -i $par_input \
    -r $par_refgene \
> $par_output