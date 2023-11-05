#!/bin/bash

set -eo pipefail

cat $par_read_1 > $par_fastq_1
cat $par_read_2 > $par_fastq_2
