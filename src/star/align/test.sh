#!/bin/bash

echo ">>> Testing $meta_functionality_name"

gunzip $meta_resources_dir/genes.gtf.gz

samtools faidx $meta_resources_dir/genome.fasta
NUM_BASES=`gawk '{sum = sum + \$2}END{if ((log(sum)/log(2))/2 - 1 > 14) {printf "%.0f", 14} else {printf "%.0f", (log(sum)/log(2))/2 - 1}}' $meta_resources_dir/genome.fasta.fai`

STAR \
    --runMode genomeGenerate \
    --genomeDir STAR_index \
    --genomeFastaFiles $meta_resources_dir/genome.fasta \
    --sjdbGTFfile $meta_resources_dir/genes.gtf \
    --genomeSAindexNbases $NUM_BASES

echo ">>> Testing for paired-end reads"

"$meta_executable" \
    --paired true \
    --input $meta_resources_dir/SRR6357070_1.fastq.gz,$meta_resources_dir/SRR6357070_2.fastq.gz \
    --star_index STAR_index \
    --gtf $meta_resources_dir/genes.gtf \
    --star_ignore_sjdbgtf false \
    --extra_star_align_args  "--readFilesCommand gunzip -c --quantMode TranscriptomeSAM --twopassMode Basic --outSAMtype BAM Unsorted --runRNGseed 0 --outFilterMultimapNmax 20 --alignSJDBoverhangMin 1 --outSAMattributes NH HI AS NM MD --quantTranscriptomeBan Singleend --outSAMstrandField intronMotif" \
    --output STAR_pe \
    --star_align_bam star_aligned_pe.genome.bam \
    --star_align_bam_transcriptome star_aligned_pe.transcriptome.bam \
    --log_final star_aligned_pe.log.final.out

echo ">> Checking if the correct files are present"
[ ! -d "STAR_pe" ] && echo "Directory 'STAR' does not exist!" && exit 1
[ -z "$(ls -A 'STAR_pe')" ] && echo "Directory 'STAR' is empty!" && exit 1
[ ! -f "star_aligned_pe.genome.bam" ] && echo "File 'star_aligned_pe.genome.bam' does not exist!" && exit 1
[ ! -s "star_aligned_pe.genome.bam" ] && echo "File 'star_aligned_pe.genome.bam' is empty!" && exit 1
[ ! -f "star_aligned_pe.transcriptome.bam" ] && echo "File 'star_aligned_pe.transcriptome.bam' does not exist!" && exit 1
[ ! -s "star_aligned_pe.transcriptome.bam" ] && echo "File 'star_aligned_pe.transcriptome.bam' is empty!" && exit 1
[ ! -f "star_aligned_pe.log.final.out" ] && echo "File 'star_aligned_pe.log.final.out' does not exist!" && exit 1
[ ! -s "star_aligned_pe.log.final.out" ] && echo "File 'star_aligned_pe.log.final.out' is empty!" && exit 1

echo ">>> Testing for single-end reads"
"$meta_executable" \
    --paired true \
    --input $meta_resources_dir/SRR6357070_1.fastq.gz \
    --star_index STAR_index \
    --gtf $meta_resources_dir/genes.gtf \
    --star_ignore_sjdbgtf false \
    --extra_star_align_args  "--readFilesCommand gunzip -c --quantMode TranscriptomeSAM --twopassMode Basic --outSAMtype BAM Unsorted --runRNGseed 0 --outFilterMultimapNmax 20 --alignSJDBoverhangMin 1 --outSAMattributes NH HI AS NM MD --quantTranscriptomeBan Singleend --outSAMstrandField intronMotif" \
    --output STAR_se \
    --star_align_bam star_aligned_se.genome.bam \
    --star_align_bam_transcriptome star_aligned_se.transcriptome.bam \
    --log_final star_aligned_se.log.final.out
    
echo ">> Checking if the correct files are present"
[ ! -d "STAR_se" ] && echo "Directory 'STAR_se' does not exist!" && exit 1
[ -z "$(ls -A 'STAR_se')" ] && echo "Directory 'STAR_se' is empty!" && exit 1
[ ! -f "star_aligned_se.genome.bam" ] && echo "File 'star_aligned_se.genome.bam' does not exist!" && exit 1
[ ! -s "star_aligned_se.genome.bam" ] && echo "File 'star_aligned_se.genome.bam' is empty!" && exit 1
[ ! -f "star_aligned_se.transcriptome.bam" ] && echo "File 'star_aligned_se.transcriptome.bam' does not exist!" && exit 1
[ ! -s "star_aligned_se.transcriptome.bam" ] && echo "File 'star_aligned_se.transcriptome.bam' is empty!" && exit 1
[ ! -f "star_aligned_se.log.final.out" ] && echo "File 'star_aligned_se.log.final.out' does not exist!" && exit 1
[ ! -s "star_aligned_se.log.final.out" ] && echo "File 'star_aligned_se.log.final.out' is empty!" && exit 1

echo ">>> Test finished successfully"
exit 0
