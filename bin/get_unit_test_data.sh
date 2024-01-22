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

# MultiQC test resources
wget https://multiqc.info/examples/rna-seq/data.zip

# RSeQC
wget https://data.cyverse.org/dav-anon/iplant/home/liguow/RSeQC/Pairend_StrandSpecific_51mer_Human_hg19.bam
wget -O hg19_RefSeq.bed.gz https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg19_RefSeq.bed.gz/download
wget https://data.cyverse.org/dav-anon/iplant/home/liguow/RSeQC/Pairend_StrandSpecific_51mer_Human_hg19.bam.bai

### Resources from https://github.com/snakemake/snakemake-wrappers/tree/master/bio
# genomecov
wget https://github.com/snakemake/snakemake-wrappers/raw/master/bio/bedtools/genomecov/test/bam_input/a.sorted.bam

# picard markduplicates
wget https://github.com/snakemake/snakemake-wrappers/raw/master/bio/picard/markduplicates/test/mapped/a.bam
wget https://github.com/snakemake/snakemake-wrappers/raw/master/bio/picard/markduplicates/test/mapped/a.bam.bai
wget https://github.com/snakemake/snakemake-wrappers/raw/master/bio/picard/markduplicates/test/ref/genome.fasta
wget https://github.com/snakemake/snakemake-wrappers/raw/master/bio/picard/markduplicates/test/ref/genome.fasta.fai

# ucsc
wget https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/genome/genome.sizes
wget https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/illumina/bedgraph/test.bedgraph

# DESeq2
wget https://github.com/snakemake/snakemake-wrappers/raw/master/bio/deseq2/deseqdataset/test/dataset/counts.tsv

# dupRadar
wget https://github.com/ssayols/dupRadar/raw/master/inst/extdata/genes.gtf
wget https://github.com/ssayols/dupRadar/raw/master/inst/extdata/wgEncodeCaltechRnaSeqGm12878R1x75dAlignsRep2V2.bam
wget https://github.com/ssayols/dupRadar/raw/master/inst/extdata/wgEncodeCaltechRnaSeqGm12878R1x75dAlignsRep2V2.bam.bai

# Qualimap
wget -O qualimap_test.bam.bai https://github.com/snakemake/snakemake-wrappers/raw/master/bio/qualimap/rnaseq/test/mapped/a.bai
wget -O qualimap_test.bam https://github.com/snakemake/snakemake-wrappers/raw/master/bio/qualimap/rnaseq/test/mapped/a.bam
wget -O qualimap_test_annot.gtf https://github.com/snakemake/snakemake-wrappers/raw/master/bio/qualimap/rnaseq/test/annotation.gtf

# preseq lc_extrap
wget https://github.com/snakemake/snakemake-wrappers/raw/master/bio/preseq/lc_extrap/test/samples/a.sorted.bed
wget https://github.com/smithlabcode/preseq/raw/master/data/SRR1106616_5M_subset.bam

cd $CURR