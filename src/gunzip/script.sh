#!/bin/bash

set -eo pipefail

if [ ${par_input##*.} == "gz" ]; then
    gunzip -c $par_input > $par_output
else
    cat $par_input > $par_output
fi