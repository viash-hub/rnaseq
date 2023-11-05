#!/usr/bin/env python3

class ExceedMaxContigSizeError(Exception):
 
    def __init__(self, value):
        self.value = value
 
    def __str__(self):
        return(repr(self.value))



## VIASH START
par = {
    "input": "testData/test_output/ref.prepare_genome.fasta_uncompressed/genome_gfp.fasta.fai"
}
## VIASH END

try:  
    max_size = 512000000
    with open(par['input']) as fai_file: 
        for line in fai_file.readlines(): 
            lspl  = line.split('\t') 
            chrom = lspl[0]
            size  = lspl[1]
            if (int(size) > max_size):
                error_string = f"""Contig longer than {max_size}bp found in reference genome!\n
                {chrom}: {size}\n
                Provide the '--bam_csi_index' parameter to use a CSI instead of BAI index.\n
                Please see: https://github.com/nf-core/rnaseq/issues/744\n"""
                raise(ExceedMaxContigSizeError(error_string))

except ExceedMaxContigSizeError as error:
    print(error.value)