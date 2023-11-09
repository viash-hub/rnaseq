#!/bin/bash

set -eo pipefail

IFS="," read -ra input <<< "$par_input"

if $star_ignore_sjdbgtf; then
    ignore_gtf='' 
else
    ignore_gtf="--sjdbGTFfile $par_gtf"
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

mkdir -p $par_output

STAR \
    --genomeDir $par_star_index \
    --readFilesIn $input \
    --runThreadN $meta_cpus \
    --outFileNamePrefix $par_output/ \
    $out_sam_type \
    $ignore_gtf \
    $attrRG \
    $par_extra_star_align_args

if [[ $par_extra_star_align_args == *"--outSAMtype BAM Unsorted SortedByCoordinate"* ]]; then
    mv ${par_output}/Aligned.out.bam ${par_output}/Aligned.unsort.out.bam 
fi

if [ -f "$par_output/Aligned.out.bam" ]; then
    cp "$par_output/Aligned.out.bam" $par_star_align_bam
elif [ -f "$par_output/Aligned.sortedByCoord.out.bam" ]; then
    cp "$par_output/Aligned.sortedByCoord.out.bam" $par_star_align_bam
fi

if [ -f "$par_output/Aligned.toTranscriptome.out.bam" ]; then
    cp "$par_output/Aligned.toTranscriptome.out.bam" $par_star_align_bam_transcriptome
fi

if [ -f "$par_output/Log.final.out" ]; then
    cp "$par_output/Log.final.out" $par_log_final
fi

if [ -f "$par_output/Unmapped.out.mate1" ]; then
    mv $par_output/Unmapped.out.mate1 $par_output/unmapped_1.fastq
    gzip $par_output/unmapped_1.fastq
fi
if [ -f "$par_output/Unmapped.out.mate2" ]; then
    mv $par_output/Unmapped.out.mate2 $par_output/unmapped_2.fastq
    gzip $par_output/unmapped_2.fastq
fi
