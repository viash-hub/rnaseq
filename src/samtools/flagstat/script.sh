#!/bin/bash

set -eo pipefail

samtools flagstat --threads $meta_cpus $par_bam > $par_output