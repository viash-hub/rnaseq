#!/bin/bash

CURR=`pwd`
DEST="testData/test/"
mkdir -p $DEST
cd $DEST

# echo "Fetching test data from https://github.com/nf-core/test-datasets/tree/rnaseq3/testdata/GSE110004..."

wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/testdata/GSE110004/SRR6357070_1.fastq.gz
wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/testdata/GSE110004/SRR6357070_2.fastq.gz

wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/testdata/GSE110004/SRR6357071_1.fastq.gz
wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/testdata/GSE110004/SRR6357071_2.fastq.gz

wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/testdata/GSE110004/SRR6357072_1.fastq.gz
wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/testdata/GSE110004/SRR6357072_2.fastq.gz

wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/testdata/GSE110004/SRR6357073_1.fastq.gz
wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/testdata/GSE110004/SRR6357074_1.fastq.gz
wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/testdata/GSE110004/SRR6357075_1.fastq.gz

wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/testdata/GSE110004/SRR6357076_1.fastq.gz
wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/testdata/GSE110004/SRR6357076_2.fastq.gz

wget https://data.cyverse.org/dav-anon/iplant/home/liguow/RSeQC/Pairend_StrandSpecific_51mer_Human_hg19.bam
wget http://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg19_rRNA.bed.gz
gunzip hg19rRNA.bed.gz

cd $CURR
