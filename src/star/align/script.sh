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
    --readFilesIn ${input[*]} \
    --runThreadN ${meta_cpus:-1} \
    --outFileNamePrefix output/ \
    $out_sam_type \
    $ignore_gtf \
    $attrRG \
    $par_extra_star_align_args

find "output/" -type f -exec cp {} "$par_output" \;

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

# Version
text="${meta_functionality_name}:
    star: $(STAR --version | sed -e "s/STAR_//g")
    samtools: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    gawk: $(echo $(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')"

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
    mv "$par_versions" "$par_updated_versions"
else
    echo "$text" > "$par_updated_versions"
fi