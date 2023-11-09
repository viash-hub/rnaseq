#!/bin/bash

set -eo pipefail

bedGraphToBigWig \
    $par_bedgraph \
    $par_sizes \
    $par_bigWig