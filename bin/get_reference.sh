#!/bin/bash

CURR=`pwd`
DEST="testData/reference/"
mkdir -p $DEST
cd $DEST

# echo "Fetching test data from https://github.com/nf-core/test-datasets/tree/rnaseq3/reference..."

wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/reference/bbsplit_fasta_list.txt
wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/reference/genes.gff.gz
wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/reference/genes.gtf.gz
wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/reference/genome.fasta
wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/reference/gfp.fa.gz
wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/reference/hisat2.tar.gz
wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/reference/rsem.tar.gz
wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/reference/salmon.tar.gz
wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/reference/transcriptome.fasta
wget https://raw.githubusercontent.com/nf-core/rnaseq/3.12.0/assets/rrna-db-defaults.txt

cd $CURR
