#!/bin/bash

CURR=`pwd`

### Get input fastq files for the minimal test

DEST_FASTQ="testData/minimal_test/input_fastq"
mkdir -p $DEST_FASTQ
cd $DEST_FASTQ

echo "Fetching FastQ files..."

# wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/testdata/GSE110004/SRR6357070_1.fastq.gz
# wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/testdata/GSE110004/SRR6357070_2.fastq.gz

# wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/testdata/GSE110004/SRR6357071_1.fastq.gz
# wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/testdata/GSE110004/SRR6357071_2.fastq.gz

# wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/testdata/GSE110004/SRR6357072_1.fastq.gz
# wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/testdata/GSE110004/SRR6357072_2.fastq.gz

# wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/testdata/GSE110004/SRR6357073_1.fastq.gz

# wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/testdata/GSE110004/SRR6357074_1.fastq.gz

# wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/testdata/GSE110004/SRR6357075_1.fastq.gz

# wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/testdata/GSE110004/SRR6357076_1.fastq.gz
# wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/testdata/GSE110004/SRR6357076_2.fastq.gz

wget https://github.com/nf-core/test-datasets/raw/rnaseq/testdata/GSE110004/SRR6357070_1.fastq.gz
wget https://github.com/nf-core/test-datasets/raw/rnaseq/testdata/GSE110004/SRR6357070_2.fastq.gz

wget https://github.com/nf-core/test-datasets/raw/rnaseq/testdata/GSE110004/SRR6357071_1.fastq.gz
wget https://github.com/nf-core/test-datasets/raw/rnaseq/testdata/GSE110004/SRR6357071_2.fastq.gz

wget https://github.com/nf-core/test-datasets/raw/rnaseq/testdata/GSE110004/SRR6357072_1.fastq.gz
wget https://github.com/nf-core/test-datasets/raw/rnaseq/testdata/GSE110004/SRR6357072_2.fastq.gz

wget https://github.com/nf-core/test-datasets/raw/rnaseq/testdata/GSE110004/SRR6357073_1.fastq.gz

wget https://github.com/nf-core/test-datasets/raw/rnaseq/testdata/GSE110004/SRR6357074_1.fastq.gz

wget https://github.com/nf-core/test-datasets/raw/rnaseq/testdata/GSE110004/SRR6357075_1.fastq.gz

wget https://github.com/nf-core/test-datasets/raw/rnaseq/testdata/GSE110004/SRR6357076_1.fastq.gz
wget https://github.com/nf-core/test-datasets/raw/rnaseq/testdata/GSE110004/SRR6357076_2.fastq.gz

cd $CURR

### Get reference files for the minimal test

DEST_REF="testData/minimal_test/reference"
mkdir -p $DEST_REF
cd $DEST_REF

echo "Fetching reference data..."

# wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/reference/bbsplit_fasta_list.txt
# wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/reference/genes.gff.gz
# wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/reference/genes.gtf.gz
# wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/reference/genome.fasta
# wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/reference/gfp.fa.gz
# wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/reference/hisat2.tar.gz
# wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/reference/rsem.tar.gz
# wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/reference/salmon.tar.gz
# wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq3/reference/transcriptome.fasta

wget https://raw.githubusercontent.com/nf-core/rnaseq/3.12.0/assets/rrna-db-defaults.txt

wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/reference/genome.fasta
wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/reference/genes.gtf.gz
wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/reference/genes.gff.gz
wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/reference/transcriptome.fasta
wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/reference/gfp.fa.gz

wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/reference/bbsplit_fasta_list.txt
# wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/reference/hisat2.tar.gz
wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/reference/salmon.tar.gz
wget https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/reference/rsem.tar.gz

cd $CURR

NEWDEST1_REF="$CURR/testData/minimal_test/reference/rRNA"
mkdir -p $NEWDEST1_REF
cd $NEWDEST1_REF
for LINE in `cat ../rrna-db-defaults.txt`
do
    wget $LINE
done
cd $CURR
find $NEWDEST1_REF -type f > $DEST_REF/rrna-db-defaults.txt

NEWDEST2_REF="$CURR/testData/minimal_test/reference/bbsplit_fasta"
mkdir -p $NEWDEST2_REF
while IFS=, read -r -a line; do
    url="${line[1]}"
    name="$NEWDEST2_REF/${line[0]}.fa"
    wget $url -O "$name"
    line+=("$name")
    IFS=','
    echo "${line[*]}" >> "$NEWDEST2_REF/tmp.txt"
done < "$DEST_REF/bbsplit_fasta_list.txt"
cut -d',' -f1,3 "$NEWDEST2_REF/tmp.txt" > "$DEST_REF/bbsplit_fasta_list.txt"
rm "$NEWDEST2_REF/tmp.txt"
