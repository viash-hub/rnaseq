#!/bin/bash

set -eo pipefail

preseq lc_extrap \
    $par_extra_preseq_args \
    ${par_paired:+-pe} \
    -output $par_output \
    $par_bam