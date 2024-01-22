#!/bin/bash

echo ">>> Testing $meta_functionality_name"

"$meta_executable" \
    --strandedness reverse \
    --bam $meta_resources_dir/wgEncodeCaltechRnaSeqGm12878R1x75dAlignsRep2V2.bam \
    --annotation_gtf $meta_resources_dir/genes.gtf \
    --extra_stringtie_args "-v" \
    --transcript_gtf test.transcripts.gtf \
    --coverage_gtf test.coverage.gtf \
    --abundance  test.abundance.txt \
    --ballgown test.ballgown 
    
echo ">> Checking if the correct files are present"
[ ! -d "test.ballgown" ] && echo "Directory 'test.ballgown' does not exist!" && exit 1
[ -z "$(ls -A 'test.ballgown')" ] && echo "Directory 'test.ballgown' is empty!" && exit 1
[ ! -f "test.transcripts.gtf" ] && echo "File 'test.transcripts.gtf' does not exist!" && exit 1
[ ! -s "test.transcripts.gtf" ] && echo "File 'test.transcripts.gtf' is empty!" && exit 1
[ ! -f "test.coverage.gtf" ] && echo "File 'test.coverage.gtf' does not exist!" && exit 1
[ ! -s "test.coverage.gtf" ] && echo "File 'test.coverage.gtf' is empty!" && exit 1
[ ! -f "test.abundance.txt" ] && echo "File 'test.abundance.txt' does not exist!" && exit 1
[ ! -s "test.abundance.txt" ] && echo "File 'test.abundance.txt' is empty!" && exit 1

echo ">>> Test finished successfully"
exit 0
