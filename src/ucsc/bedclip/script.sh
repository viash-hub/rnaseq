#!/bin/bash

set -eo pipefail

bedClip \
    $par_input_bedgraph \
    $par_sizes \
    $par_output_bedGraph