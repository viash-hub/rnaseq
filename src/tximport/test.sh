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

cat > tx2gene.tsv << HERE
ENSSASG00005000002      ENSSASG00005000002      ORF1ab
ENSSASG00005000003      ENSSASG00005000003      ORF1ab
ENSSASG00005000004      ENSSASG00005000004      S
HERE

echo "> Run test"
"$meta_executable" \
    --quant_results "sample1_quant_results.sf;sample2_quant_results.sf" \
    --tx2gene_tsv "tx2gene.tsv" \
    --quant_type "salmon" \
    # --tpm_gene gene_tpm.tsv \
    # --counts_gene gene_counts.tsv \
    # --counts_gene_length_scaled gene_counts_length_scaled.tsv \
    # --counts_gene_scaled gene_counts_scaled.tsv \
    # --lengths_gene gene_length.tsv \
    # --tpm_transcript transcript_tpm.tsv \
    # --counts_transcript transcript_counts.tsv \
    # --lengths_transcript transcript_length.tsv

echo "> Check results"
[ ! -f "merged.gene_tpm.tsv" ] && echo "merged.gene_tpm.tsv was not created" && exit 1
[ ! -s "merged.gene_tpm.tsv" ] && echo "merged.gene_tpm.tsv is empty" && exit 1
[ ! -f "merged.gene_counts.tsv" ] && echo "merged.gene_counts.tsv was not created" && exit 1
[ ! -s "merged.gene_counts.tsv" ] && echo "merged.gene_counts.tsv is empty" && exit 1
[ ! -f "merged.gene_counts_length_scaled.tsv" ] && echo "merged.gene_counts_length_scaled.tsv was not created" && exit 1
[ ! -s "merged.gene_counts_length_scaled.tsv" ] && echo "merged.gene_counts_length_scaled.tsv is empty" && exit 1
[ ! -f "merged.gene_counts_scaled.tsv" ] && echo "merged.gene_counts_scaled.tsv was not created" && exit 1
[ ! -s "merged.gene_counts_scaled.tsv" ] && echo "merged.gene_counts_scaled.tsv is empty" && exit 1
[ ! -f "merged.gene_lengths.tsv" ] && echo "merged.gene_lengths.tsv was not created" && exit 1
[ ! -s "merged.gene_lengths.tsv" ] && echo "merged.gene_lengths.tsv is empty" && exit 1
[ ! -f "merged.transcript_tpm.tsv" ] && echo "merged.transcript_tpm.tsv was not created" && exit 1
[ ! -s "merged.transcript_tpm.tsv" ] && echo "merged.transcript_tpm.tsv is empty" && exit 1
[ ! -f "merged.transcript_counts.tsv" ] && echo "merged.transcript_counts.tsv was not created" && exit 1
[ ! -s "merged.transcript_counts.tsv" ] && echo "merged.transcript_counts.tsv is empty" && exit 1
[ ! -f "merged.transcript_lengths.tsv" ] && echo "merged.transcript_lengths.tsv was not created" && exit 1
[ ! -s "merged.transcript_lengths.tsv" ] && echo "merged.transcript_lengths.tsv is empty" && exit 1

exit 0