#!/bin/bash

set -e 

echo "> Prepare test data"

cat > "test.featurecounts.txt" << HERE
# Program:featureCounts v2.0.1; Command:"featureCounts" "-t" "CDS" "-T" "2" "-a" "genome.gtf" "-s" "1" "-o" "test.featureCounts.txt" "test.single_end.bam"
Geneid	Chr	Start	End	Strand	Length	test.single_end.bam
orf1ab	MT192765.1;MT192765.1	259;13461	13461;21545	+;+	21287	38
S	MT192765.1	21556	25374	+	3819	4
ORF3a	MT192765.1	25386	26210	+	825	0
E	MT192765.1	26238	26462	+	225	1
M	MT192765.1	26516	27181	+	666	1
ORF6	MT192765.1	27195	27377	+	183	0
ORF7a	MT192765.1	27387	27749	+	363	0
ORF7b	MT192765.1	27749	27877	+	129	0
ORF8	MT192765.1	27887	28249	+	363	0
N	MT192765.1	28267	29523	+	1257	2
ORF10	MT192765.1	29551	29664	+	114	0
HERE

echo "> Run test"
"$meta_executable" \
    --biocounts "test.featurecounts.txt" \
    --id "test" \
    --featurecounts_multiqc "test.biotype_counts_mqc.tsv" \
    --featurecounts_rrna_multiqc "test.biotype_counts_rrna_mqc.tsv"

echo "> Check results"
[ ! -f "test.biotype_counts_mqc.tsv" ] && echo "test.biotype_counts_mqc.tsv was not created" && exit 1
[ ! -s "test.biotype_counts_mqc.tsv" ] && echo "test.biotype_counts_mqc.tsv is empty" && exit 1
[ ! -f "test.biotype_counts_rrna_mqc.tsv" ] && echo "test.biotype_counts_rrna_mqc.tsv was not created" && exit 1
[ ! -s "test.biotype_counts_rrna_mqc.tsv" ] && echo "test.biotype_counts_rrna_mqc.tsv is empty" && exit 1

exit 0
