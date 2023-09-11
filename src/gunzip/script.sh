#!/bin/bash

set -eo pipefail

if [ "${par_input##*.}" == "gz" ]; then
    echo "${par_input##*.}"
    gunzip -f "$par_input" > "$par_output"
else
    cat $par_input > $par_output
fi