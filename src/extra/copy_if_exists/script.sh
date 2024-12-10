#!/bin/bash

set -eo pipefail

mkdir -p $par_output

# This file is checked by the Nextflow module wrapper
cp $par_required_file "$par_output"

# If the variable is empty, we use the default one (registered as a resource)
if [ -z $par_optional_file ]; then
  echo "No optional_file provided, using the default"
  cp $meta_resources_dir/optional_file.txt "$par_output"
else
  echo "Optional file provided"
  if [ -f $par_optional_file ]; then
    cp $par_optional_file "$par_output"
  else
    # Unreachable: the Viash-generated module checks this
    echo "Optional file does not exist"
    exit 1
  fi
fi

echo "Done"
