#!/bin/bash

set -eo pipefail

samtools sort -@ $meta_cpus -o $par_output $par_input