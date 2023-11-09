#!/bin/bash

set -eo pipefail

samtools idxstats --threads $meta_cpus $par_bam > $par_output