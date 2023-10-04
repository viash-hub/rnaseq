#!/bin/bash

set -eo pipefail

## VIASH START
# par_fastq_1="testData/paired_end_test/SRR6357070.sortmerna.output/non_rRNA_reads_1.fq.gz" 
# par_fastq_2="testData/paired_end_test/SRR6357070.sortmerna.output/non_rRNA_reads_2.fq.gz"
par_paired=true, 
input=("testData/test/SRR6357070_1.fastq.gz" "testData/test/SRR6357070_2.fastq.gz")
par_gtf="testData/test_output/ref.gtf_gene_filter.filtered_gtf"
par_star_index="testData/test_output/ref.star_index_uncompressed.star_index"
par_star_ignore_sjdbgtf=false
par_seq_platform=''
par_seq_center=''
par_extra_star_align_args='--readFilesCommand gunzip -c'
meta_cpus=2
par_output="test"

## VIASH END

# IFS="," read -ra input <<< "$par_input"

# def reads1 = [], reads2 = []
# meta.single_end ? [reads].flatten().each{reads1 << it} : reads.eachWithIndex{ v, ix -> ( ix & 1 ? reads2 : reads1) << v }

if $star_ignore_sjdbgtf; then
    ignore_gtf='' 
else
    ignore_gtf="--sjdbGTFfile $gtf"
fi

if $par_seq_platform; then
    seq_platform="PL:$par_seq_platform" 
else
    seq_platform=''
fi

if $par_seq_center; then
    seq_center="CN:$par_seq_center" 
else
    seq_center=''
fi

if [[ $par_extra_star_align_args == *"--outSAMattrRGline"* ]]; then
    attrRG="" 
else
    attrRG="--outSAMattrRGline ID:$par_output $seq_center SM:$par_output $seq_platform"
fi

if [[ $par_extra_star_align_args == *"--outSAMtype"* ]]; then
    out_sam_type="" 
else
    out_sam_type="--outSAMtype BAM Unsorted"
fi

if [[ $par_extra_star_align_args == *"--outSAMtype BAM Unsorted SortedByCoordinate"* ]]; then
    mv_unsorted_bam="mv ${par_output}/Aligned.out.bam ${par_output}/Aligned.unsort.out.bam" 
else
    mv_unsorted_bam=""
fi

mkdir -p $par_output
STAR --genomeDir $par_star_index --readFilesIn ${input[0]} ${input[1]} --runThreadN $meta_cpus --outFileNamePrefix $par_output/ $out_sam_type $ignore_gtf $attrRG $par_extra_star_align_args

$mv_unsorted_bam

if [ -f "$par_output/Unmapped.out.mate1" ]; then
    mv $par_output/Unmapped.out.mate1 $par_output/unmapped_1.fastq
    gzip $par_output/unmapped_1.fastq
fi
if [ -f "$par_output/Unmapped.out.mate2" ]; then
    mv $par_output/Unmapped.out.mate2 $par_output/unmapped_2.fastq
    gzip $par_output/unmapped_2.fastq
fi
