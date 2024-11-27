#!/bin/bash 

set -e

echo "> Prepare test data"

cat > "sample1_quant_results.sf" << HERE
Name	Length	EffectiveLength	TPM	NumReads
ENSSASG00005000004	3822	3572	15216.8	753
ENSSASG00005000003	13218	12968	1502.34	269.9
ENSSASG00005000002	21290	21040	23916.3	6971.1
HERE

cat > "sample2_quant_results.sf" << HERE
Name	Length	EffectiveLength	TPM	NumReads
ENSSASG00005000004	3822	3572	23713.5	703
ENSSASG00005000003	13218	12968	14280	1536.92
ENSSASG00005000002	21290	21040	37447.4	6539.08
HERE

cat > "genes.gtf" << HERE
MN908947.3	ensembl	transcript	266	21555	.	+	.	gene_id "ENSSASG00005000002"; gene_version "1"; transcript_id "ENSSAST00005000002"; transcript_version "1"; gene_name "ORF1ab"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "ORF1ab"; transcript_source "ensembl"; transcript_biotype "protein_coding";
MN908947.3	ensembl	exon	266	21555	.	+	.	gene_id "ENSSASG00005000002"; gene_version "1"; transcript_id "ENSSAST00005000002"; transcript_version "1"; exon_number "1"; gene_name "ORF1ab"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "ORF1ab"; transcript_source "ensembl"; transcript_biotype "protein_coding"; exon_id "ENSSASE00005000002"; exon_version "1";
MN908947.3	ensembl	CDS	266	21552	.	+	0	gene_id "ENSSASG00005000002"; gene_version "1"; transcript_id "ENSSAST00005000002"; transcript_version "1"; exon_number "1"; gene_name "ORF1ab"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "ORF1ab"; transcript_source "ensembl"; transcript_biotype "protein_coding"; protein_id "ENSSASP00005000002"; protein_version "1";
MN908947.3	ensembl	start_codon	266	268	.	+	0	gene_id "ENSSASG00005000002"; gene_version "1"; transcript_id "ENSSAST00005000002"; transcript_version "1"; exon_number "1"; gene_name "ORF1ab"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "ORF1ab"; transcript_source "ensembl"; transcript_biotype "protein_coding";
MN908947.3	ensembl	stop_codon	21553	21555	.	+	0	gene_id "ENSSASG00005000002"; gene_version "1"; transcript_id "ENSSAST00005000002"; transcript_version "1"; exon_number "1"; gene_name "ORF1ab"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "ORF1ab"; transcript_source "ensembl"; transcript_biotype "protein_coding";
MN908947.3	ensembl	transcript	266	13483	.	+	.	gene_id "ENSSASG00005000003"; gene_version "1"; transcript_id "ENSSAST00005000003"; transcript_version "1"; gene_name "ORF1ab"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "ORF1a"; transcript_source "ensembl"; transcript_biotype "protein_coding";
MN908947.3	ensembl	exon	266	13483	.	+	.	gene_id "ENSSASG00005000003"; gene_version "1"; transcript_id "ENSSAST00005000003"; transcript_version "1"; exon_number "1"; gene_name "ORF1ab"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "ORF1a"; transcript_source "ensembl"; transcript_biotype "protein_coding"; exon_id "ENSSASE00005000003"; exon_version "1";
MN908947.3	ensembl	CDS	266	13480	.	+	0	gene_id "ENSSASG00005000003"; gene_version "1"; transcript_id "ENSSAST00005000003"; transcript_version "1"; exon_number "1"; gene_name "ORF1ab"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "ORF1a"; transcript_source "ensembl"; transcript_biotype "protein_coding"; protein_id "ENSSASP00005000003"; protein_version "1";
MN908947.3	ensembl	start_codon	266	268	.	+	0	gene_id "ENSSASG00005000003"; gene_version "1"; transcript_id "ENSSAST00005000003"; transcript_version "1"; exon_number "1"; gene_name "ORF1ab"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "ORF1a"; transcript_source "ensembl"; transcript_biotype "protein_coding";
MN908947.3	ensembl	stop_codon	13481	13483	.	+	0	gene_id "ENSSASG00005000003"; gene_version "1"; transcript_id "ENSSAST00005000003"; transcript_version "1"; exon_number "1"; gene_name "ORF1ab"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "ORF1a"; transcript_source "ensembl"; transcript_biotype "protein_coding";
MN908947.3	ensembl	transcript	21563	25384	.	+	.	gene_id "ENSSASG00005000004"; gene_version "1"; transcript_id "ENSSAST00005000004"; transcript_version "1"; gene_name "S"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "S"; transcript_source "ensembl"; transcript_biotype "protein_coding";
MN908947.3	ensembl	exon	21563	25384	.	+	.	gene_id "ENSSASG00005000004"; gene_version "1"; transcript_id "ENSSAST00005000004"; transcript_version "1"; exon_number "1"; gene_name "S"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "S"; transcript_source "ensembl"; transcript_biotype "protein_coding"; exon_id "ENSSASE00005000004"; exon_version "1";
MN908947.3	ensembl	CDS	21563	25381	.	+	0	gene_id "ENSSASG00005000004"; gene_version "1"; transcript_id "ENSSAST00005000004"; transcript_version "1"; exon_number "1"; gene_name "S"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "S"; transcript_source "ensembl"; transcript_biotype "protein_coding"; protein_id "ENSSASP00005000004"; protein_version "1";
MN908947.3	ensembl	start_codon	21563	21565	.	+	0	gene_id "ENSSASG00005000004"; gene_version "1"; transcript_id "ENSSAST00005000004"; transcript_version "1"; exon_number "1"; gene_name "S"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "S"; transcript_source "ensembl"; transcript_biotype "protein_coding";
MN908947.3	ensembl	stop_codon	25382	25384	.	+	0	gene_id "ENSSASG00005000004"; gene_version "1"; transcript_id "ENSSAST00005000004"; transcript_version "1"; exon_number "1"; gene_name "S"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "S"; transcript_source "ensembl"; transcript_biotype "protein_coding";
HERE

echo "> Run test"
"$meta_executable" \
    --quant_results "sample1_quant_results.sf;sample2_quant_results.sf" \
    --gtf "genes.gtf" \
    --quant_type "salmon" \
    --gtf_extra_attributes "gene_name" \
    --gtf_group_features "gene_id" \
    --tsv "tx2gene.tsv" 

echo "> Check results"
[ ! -f "tx2gene.tsv" ] && echo "tx2gene.tsv was not created" && exit 1
[ ! -s "tx2gene.tsv" ] && echo "tx2gene.tsv is empty" && exit 1

exit 0
