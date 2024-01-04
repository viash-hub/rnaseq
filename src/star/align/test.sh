#!/bin/bash

echo ">>> Testing $meta_functionality_name"

gunzip $meta_resources_dir/genes.gtf.gz

echo ">>> Testing for paired-end reads"

"$meta_executable" \
    --paired true \
    --input $meta_resources_dir/SRR6357070_1.fastq.gz,$meta_resources_dir/SRR6357070_2.fastq.gz \
    --gtf $meta_resources_dir/genes.gtf \
    --star_ignore_sjdbgtf false \
    --extra_star_align_args  "--readFilesCommand gunzip -c --quantMode TranscriptomeSAM --twopassMode Basic --outSAMtype BAM Unsorted --runRNGseed 0 --outFilterMultimapNmax 20 --alignSJDBoverhangMin 1 --outSAMattributes NH HI AS NM MD --quantTranscriptomeBan Singleend --outSAMstrandField intronMotif" \
    --output  STAR \
    --star_align_bam star_aligned.genome.bam \
    --star_align_bam_transcriptome star_aligned.transcriptome.bam \
    --log_final star_aligned.log.final.out

echo ">> Checking if the correct files are present"
[ ! -d STAR ] && echo "Directory 'STAR' does not exist!" && exit 1
[ ! -f star_aligned.genome.bam ] && echo "File 'star_aligned.genome.bam' does not exist!" && exit 1
[ ! -f star_aligned.transcriptome.bam ] && echo "File 'star_aligned.transcriptome.bam' does not exist!" && exit 1
[ ! -f star_aligned.log.final.out ] && echo "File 'star_aligned.log.final.out' does not exist!" && exit 1

echo ">>> Testing for single-end reads"
"$meta_executable" \
    --paired true \
    --input $meta_resources_dir/SRR6357070_1.fastq.gz \
    --gtf $meta_resources_dir/genes.gtf \
    --star_ignore_sjdbgtf false \
    --extra_star_align_args  "--readFilesCommand gunzip -c --quantMode TranscriptomeSAM --twopassMode Basic --outSAMtype BAM Unsorted --runRNGseed 0 --outFilterMultimapNmax 20 --alignSJDBoverhangMin 1 --outSAMattributes NH HI AS NM MD --quantTranscriptomeBan Singleend --outSAMstrandField intronMotif" \
    --output  STAR \
    --star_align_bam star_aligned.genome.bam \
    --star_align_bam_transcriptome star_aligned.transcriptome.bam \
    --log_final star_aligned.log.final.out
    
echo ">> Checking if the correct files are present"
[ ! -d STAR ] && echo "Directory 'STAR' does not exist!" && exit 1
[ ! -f star_aligned.genome.bam ] && echo "File 'star_aligned.genome.bam' does not exist!" && exit 1
[ ! -f star_aligned.transcriptome.bam ] && echo "File 'star_aligned.transcriptome.bam' does not exist!" && exit 1
[ ! -f star_aligned.log.final.out ] && echo "File 'star_aligned.log.final.out' does not exist!" && exit 1

echo ">>> Test finished successfully"
exit 0
