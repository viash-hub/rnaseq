#!/bin/bash

viash ns build -q prepare_reference --setup cb

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

nextflow \
  run . \
  -main-script src/prepare_reference/test.nf \
  -entry test_wf \
  -profile docker,no_publish,local \
  -c src/configs/labels_ci.config \
  -c src/configs/integration_tests.config \

