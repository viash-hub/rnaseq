#!/bin/bash

set -eo pipefail

samtools faidx $par_fasta cut -f 1,2 $par_fai > $par_sizes