#!/bin/bash

CURR=`pwd`
DEST="testData/unit_test_resources"
mkdir -p $DEST
cd $DEST

echo "Fetching unit test resources..."

## UMI_TOOLS
# extract
wget https://github.com/CGATOxford/UMI-tools/raw/master/tests/slim.fastq.gz
wget https://github.com/CGATOxford/UMI-tools/raw/master/tests/scrb_seq_fastq.1.gz
wget https://github.com/CGATOxford/UMI-tools/raw/master/tests/scrb_seq_fastq.2.gz
# dedup
wget https://github.com/CGATOxford/UMI-tools/raw/master/tests/chr19.bam
wget https://github.com/CGATOxford/UMI-tools/raw/master/tests/chr19.bam.bai

# MultiQC
wget https://multiqc.info/examples/rna-seq/data.zip

# dupRadar
wget https://github.com/ssayols/dupRadar/raw/master/inst/extdata/genes.gtf
wget https://github.com/ssayols/dupRadar/raw/master/inst/extdata/wgEncodeCaltechRnaSeqGm12878R1x75dAlignsRep2V2.bam
wget https://github.com/ssayols/dupRadar/raw/master/inst/extdata/wgEncodeCaltechRnaSeqGm12878R1x75dAlignsRep2V2.bam.bai


### Resources from https://github.com/snakemake/snakemake-wrappers/tree/master/bio
# DESeq2
wget https://github.com/snakemake/snakemake-wrappers/raw/master/bio/deseq2/deseqdataset/test/dataset/counts.tsv

# preseq lc_extrap
wget https://github.com/snakemake/snakemake-wrappers/raw/master/bio/preseq/lc_extrap/test/samples/a.sorted.bed
wget https://github.com/smithlabcode/preseq/raw/master/data/SRR1106616_5M_subset.bam


### nf-core test datasets
# sarscov2
mkdir -p sarscov2
wget -O sarscov2/genome.sizes https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/genome/genome.sizes
wget -O sarscov2/test.bedgraph https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/illumina/bedgraph/test.bedgraph
wget -O sarscov2/genome.fasta https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/genome/genome.fasta
wget -O sarscov2/genome.fasta.fai https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/genome/genome.fasta.fai
wget -O sarscov2/test.paired_end.sorted.bam https://github.com/nf-core/test-datasets/raw/modules/data/genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam
wget -O sarscov2/test.paired_end.sorted.bam.bai https://github.com/nf-core/test-datasets/raw/modules/data/genomics/sarscov2/illumina/bam/test.paired_end.sorted.bam.bai
wget -O sarscov2/test.bed https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/genome/bed/test.bed
wget -O sarscov2/test.bed12 https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/genome/bed/test.bed12
wget -O sarscov2/genome.gtf https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/genome/genome.gtf

cd $CURR