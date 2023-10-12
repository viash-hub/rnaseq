#!/bin/bash

CURR=`pwd`
DEST="testData/reference"
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

NEWDEST1="testData/reference/rRNA"
mkdir -p $NEWDEST1
cd $NEWDEST1
for LINE in `cat ../rrna-db-defaults.txt`
do
    wget $LINE
done
cd $CURR
find $NEWDEST1 -type f > $DEST/rrna-db-defaults.txt

NEWDEST2="testData/reference/bbsplit_fasta"
mkdir -p $NEWDEST2
while IFS=, read -r -a line; do
    url="${line[1]}"
    name="$NEWDEST2/${line[0]}.fa"
    wget $url -O "$name"
    line+=("$name")
    IFS=','
    echo "${line[*]}" >> "$NEWDEST2/tmp.txt"
done < "$DEST/bbsplit_fasta_list.txt"
cut -d',' -f1,3 "$NEWDEST2/tmp.txt" > "$DEST/bbsplit_fasta_list.txt"
rm "$NEWDEST2/tmp.txt"
