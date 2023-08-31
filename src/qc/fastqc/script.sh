#!/bin/bash

set -eo pipefail

mkdir -p $par_output 

fastqc $par_input -o $par_output
