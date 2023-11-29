#!/bin/bash

## VIASH START
meta_resources_dir="..."
meta_executable="..."
## VIASH END

echo ">>> Testing with paired-end reads"

"$meta_executable" \
  --id mysample_id \
  --paired true \
  --input "$meta_resources_dir/some_fastq/input_r1.fastq,$meta_resources_dir/some_fastq/input_r2.fastq" \
  ... other params ... \
  --bbsplit_index foo \
  --filtered_output bar

# check whether output exists
[ ! -d foo ] && "Directory 'foo' does not exist!" && exit 1
[ ! -d bar ] && "Directory 'bar' does not exist!" && exit 1

echo ">>> Testing with single-end reads"

"$meta_executable" \
  --id mysample_id \
  --paired true \
  --input "$meta_resources_dir/some_fastq/input_r1.fastq,$meta_resources_dir/some_fastq/input_r2.fastq" \
  ... other params ... \
  --bbsplit_index foo \
  --filtered_output bar

# check whether output exists
[ ! -d foo ] && "Directory 'foo' does not exist!" && exit 1
[ ! -d bar ] && "Directory 'bar' does not exist!" && exit 1

# TODO: check contents of output

echo "All tests succeeded!"
exit 0