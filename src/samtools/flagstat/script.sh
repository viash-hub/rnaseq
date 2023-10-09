#!/bin/bash

set -eo pipefail

## VIASH START
meta_cpus=2
## VIASH END

samtools flagstat --threads $meta_cpus $par_input/*.bam > $par_output