#!/bin/bash

set -eo pipefail

par_output="${par_input%%.gz}"
gunzip -f $par_input 