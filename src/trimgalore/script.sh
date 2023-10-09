#!/bin/bash

set -eo pipefail

mkdir -p $par_output

read1=$(find $par_input/ -iname *_1.f*q*)
read2=$(find $par_input/ -iname *_2.f*q*)

if [ "$par_paired" == "true" ]; then
    trim_galore $par_extra_trimgalore_args --paired --gzip $read1 $read2 -o $par_output
else
    trim_galore $par_extra_trimgalore_args --gzip $read1 -o $par_output
fi 