#!/bin/bash

set -eo pipefail 

perl "$meta_resources_dir/gtf2bed" $par_gtf > $par_bed_output