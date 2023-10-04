#!/bin/bash

set -eo pipefail

## VIASH START
par_input="test_index"
par_output="test_flagstat"
meta_cpus=2
## VIASH END

samtools flagstat --threads $meta_cpus $reference $par_input/*.bam > $par_output