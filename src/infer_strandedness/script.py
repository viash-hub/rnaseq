#!/usr/bin/env python3

import json

## VIASH START
par = {
    "json": "testData/test_output/reference_genome.fasta.fai", 
    "strandedness": ""
}
## VIASH END

with open(par['json']) as json_file:
    data = json.load(json_file)

lib_type = data['library_types'][0]
par['strandedness'] = 'reverse'
if lib_type:
    if lib_type in ['U', 'IU']:
        par['strandedness'] = 'unstranded'
    elif lib_type in ['SF', 'ISF']:
        par['strandedness'] = 'forward'
    elif lib_type in ['SR', 'ISR']:
        par['strandedness'] = 'reverse'
