#!/bin/bash

set -eo pipefail 

"$meta_resources_dir/gtf2bed" $par_gtf > $par_bed_output