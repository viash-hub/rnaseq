#!/bin/bash

set -eo pipefail
par_input="/home/nirmayi/data_intuitive/rnaseq.vsh/testData/reference/genes.gff.gz"

par_output=(echo "${par_input%%.gz}")

gunzip -f $par_input 