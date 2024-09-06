#!/bin/bash

set -eo pipefail 

infer_experiment.py \
    -i $par_input \
    -r $par_refgene \
    -s $par_sample_size \
    -q $par_map_qual \
> $par_output
