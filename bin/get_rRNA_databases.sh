#!/bin/bash

CURR=`pwd`
DEST="testData/reference/rRNA"
mkdir -p $DEST
cd $DEST

wget https://raw.githubusercontent.com/biocore/sortmerna/v4.3.4/data/rRNA_databases/rfam-5.8s-database-id98.fasta
wget https://raw.githubusercontent.com/biocore/sortmerna/v4.3.4/data/rRNA_databases/rfam-5s-database-id98.fasta
wget https://raw.githubusercontent.com/biocore/sortmerna/v4.3.4/data/rRNA_databases/silva-arc-16s-id95.fasta
wget https://raw.githubusercontent.com/biocore/sortmerna/v4.3.4/data/rRNA_databases/silva-arc-23s-id98.fasta
wget https://raw.githubusercontent.com/biocore/sortmerna/v4.3.4/data/rRNA_databases/silva-bac-16s-id90.fasta
wget https://raw.githubusercontent.com/biocore/sortmerna/v4.3.4/data/rRNA_databases/silva-bac-23s-id98.fasta
wget https://raw.githubusercontent.com/biocore/sortmerna/v4.3.4/data/rRNA_databases/silva-euk-18s-id95.fasta
wget https://raw.githubusercontent.com/biocore/sortmerna/v4.3.4/data/rRNA_databases/silva-euk-28s-id98.fasta

cd $CURR
