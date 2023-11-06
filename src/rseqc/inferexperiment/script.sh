#!/bin/bash

set -eo pipefail 

infer_experiment.py \
    -i $par_input \
    -r $par_refgene \
    -s $par_sample_size \
> $par_output