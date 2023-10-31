#!/bin/bash

set -eo pipefail

meta_cpus=2

IFS="," read -ra input <<< "$par_input"

if $star_ignore_sjdbgtf; then
    ignore_gtf='' 
else
    ignore_gtf="--sjdbGTFfile $par_gtf/*"
fi

if [[ $par_extra_star_align_args == *"--outSAMattrRGline"* ]]; then
    attrRG="" 
else
    attrRG="--outSAMattrRGline ID:$par_output CN:$seq_center SM:$par_output PL:$seq_platform"
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
mkdir -p $par_star_align_bam
mkdir -p $par_star_align_bam_transcriptome

STAR \
    --genomeDir $par_star_index \
    --readFilesIn $input \
    --runThreadN $meta_cpus \
    --outFileNamePrefix $par_output/ \
    $out_sam_type \
    $ignore_gtf \
    $attrRG \
    $par_extra_star_align_args

cp "$par_output/Aligned.out.bam" "$par_star_align_bam/aligned.genome.bam"
cp "$par_output/Aligned.toTranscriptome.out.bam" "$par_star_align_bam_transcriptome/aligned.transcriptome.bam"
cp "$par_output/Log.final.out" "$par_log_final"

$mv_unsorted_bam

if [ -f "$par_output/Unmapped.out.mate1" ]; then
    mv $par_output/Unmapped.out.mate1 $par_output/unmapped_1.fastq
    gzip $par_output/unmapped_1.fastq
fi
if [ -f "$par_output/Unmapped.out.mate2" ]; then
    mv $par_output/Unmapped.out.mate2 $par_output/unmapped_2.fastq
    gzip $par_output/unmapped_2.fastq
fi
