#!/bin/bash

set -eo pipefail 

perl "$meta_resources_dir/gtf2bed.pl" $par_gtf > $par_bed_output
