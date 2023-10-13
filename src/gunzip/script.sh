#!/bin/bash

set -eo pipefail

filename="$(basename -- "$par_input")"
mkdir -p $par_output

if [ ${filename##*.} == "gz" ]; then
    gunzip -c $par_input > $par_output/${filename%.*}
else
    cat $par_input > $par_output/$filename
fi