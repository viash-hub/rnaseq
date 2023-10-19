#!/bin/bash

## VIASH START
# meta_resources_dir="src/umitools_prepareforquant"
# par_bam="testData/paired_end_test/SRR6357070.transcriptome_deduped_index.output"
# par_output="test_componenet"
# par_log="test.log"
## VIASH END

set -eo pipefail

filename="$(basename -- $par_bam/*.bam)"
echo $filename
"$meta_resources_dir/prepare-for-rsem.py" --stdin=$par_bam/$filename --stdout=$par_output/$filename --log=$par_log 