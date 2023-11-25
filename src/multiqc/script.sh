#!/bin/bash

set -eo pipefail

function clean_up {
    rm -rf "$multiqc_tmpdir"
}
trap clean_up EXIT

multiqc_tmpdir=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXXXX")

[ -e "$par_fail_trimming_multiqc" ] && cp -r $par_fail_trimming_multiqc $multiqc_tmpdir/
[ -e "$par_fail_mapping_multiqc" ] && cp -r $par_fail_mapping_multiqc $multiqc_tmpdir/
[ -e "$par_fail_strand_multiqc" ] && cp -r $par_fail_strand_multiqc $multiqc_tmpdir/
IFS="," read -ra fastqc_raw_multiqc <<< $par_fastqc_raw_multiqc && for file in "${fastqc_raw_multiqc[@]}"; do [ -e "$file" ] && cp -r $file $multiqc_tmpdir/; done
IFS="," read -ra fastqc_trim_multiqc <<< $par_fastqc_trim_multiqc && for file in "${fastqc_trim_multiqc[@]}"; do [ -e "$file" ] && cp -r $file $multiqc_tmpdir/; done
IFS="," read -ra fastqc_trim_multiqc <<< $par_fastqc_trim_multiqc && for file in "${trim_log_multiqc[@]}"; do [ -e "$file" ] && cp -r $file $multiqc_tmpdir/; done
IFS="," read -ra sortmerna_multiqc <<< $par_sortmerna_multiqc && for file in "${sortmerna_multiqc[@]}"; do [ -e "$file" ] && cp -r $file $multiqc_tmpdir/; done
IFS="," read -ra star_multiqc <<< $par_star_multiqc && for file in "${star_multiqc[@]}"; do [ -e "$file" ] && cp -r $file $multiqc_tmpdir/; done
IFS="," read -ra hisat2_multiqc <<< $par_hisat2_multiqc && for file in "${hisat2_multiqc[@]}"; do [ -e "$file" ] && cp -r $file $multiqc_tmpdir/; done
IFS="," read -ra rsem_multiqc <<< $par_rsem_multiqc && for file in "${rsem_multiqc[@]}"; do [ -e "$file" ] && cp -r $file $multiqc_tmpdir/; done
IFS="," read -ra salmon_multiqc <<< $par_salmon_multiqc && for file in "${salmon_multiqc[@]}"; do [ -e "$file" ] && cp -r $file $multiqc_tmpdir/; done
IFS="," read -ra samtools_stats <<< $par_samtools_stats && for file in "${samtools_stats[@]}"; do [ -e "$file" ] && cp -r $file $multiqc_tmpdir/; done
IFS="," read -ra samtools_flagstat <<< $par_samtools_flagstat && for file in "${samtools_flagstat[@]}"; do [ -e "$file" ] && cp -r $file $multiqc_tmpdir/; done
IFS="," read -ra samtools_idxstats <<< $par_samtools_idxstats && for file in "${samtools_idxstats[@]}"; do [ -e "$file" ] && cp -r $file $multiqc_tmpdir/; done
IFS="," read -ra markduplicates_multiqc <<< $par_markduplicates_multiqc && for file in "${markduplicates_multiqc[@]}"; do [ -e "$file" ] && cp -r $file $multiqc_tmpdir/; done
IFS="," read -ra featurecounts_multiqc <<< $par_featurecounts_multiqc && for file in "${featurecounts_multiqc[@]}"; do [ -e "$file" ] && cp -r $file $multiqc_tmpdir/; done
[ -e "$par_aligner_pca_multiqc" ] && cp -r $par_aligner_pca_multiqc $multiqc_tmpdir/
[ -e "$par_aligner_clustering_multiqc" ] && cp -r $par_aligner_clustering_multiqc $multiqc_tmpdir/
[ -e "$par_pseudoaligner_pca_multiqc" ] && cp -r $par_pseudoaligner_pca_multiqc $multiqc_tmpdir/
[ -e "$par_pseudoaligner_clustering_multiqc" ] && cp -r $par_pseudoaligner_clustering_multiqc $multiqc_tmpdir/
IFS="," read -ra preseq_multiqc <<< $par_preseq_multiqc && for file in "${preseq_multiqc[@]}"; do [ -e "$file" ] && cp -r $file $multiqc_tmpdir/; done
IFS="," read -ra qualimap_multiqc <<< $par_qualimap_multiqc && for file in "${qualimap_multiqc[@]}"; do [ -e "$file" ] && cp -r $file $multiqc_tmpdir/; done
IFS="," read -ra dupradar_multiqc <<< $par_dupradar_multiqc && for file in "${dupradar_multiqc[@]}"; do [ -e "$file" ] && cp -r $file $multiqc_tmpdir/; done
IFS="," read -ra bamstat_multiqc <<< $par_bamstat_multiqc && for file in "${bamstat_multiqc[@]}"; do [ -e "$file" ] && cp -r $file $multiqc_tmpdir/; done
IFS="," read -ra inferexperiment_multiqc <<< $par_inferexperiment_multiqc && for file in "${inferexperiment_multiqc[@]}"; do [ -e "$file" ] && cp -r $file $multiqc_tmpdir/; done
IFS="," read -ra innerdistance_multiqc <<< $par_innerdistance_multiqc && for file in "${innerdistance_multiqc[@]}"; do [ -e "$file" ] && cp -r $file $multiqc_tmpdir/; done
IFS="," read -ra junctionannotation_multiqc <<< $par_junctionannotation_multiqc && for file in "${junctionannotation_multiqc[@]}"; do [ -e "$file" ] && cp -r $file $multiqc_tmpdir/; done
IFS="," read -ra unctionsaturation_multiqc <<< $par_junctionsaturation_multiqc && for file in "${junctionsaturation_multiqc[@]}"; do [ -e "$file" ] && cp -r $file $multiqc_tmpdir/; done
IFS="," read -ra readdistribution_multiqc <<< $par_readdistribution_multiqc && for file in "${readdistribution_multiqc[@]}"; do [ -e "$file" ] && cp -r $file $multiqc_tmpdir/; done
IFS="," read -ra readduplication_multiqc <<< $par_readduplication_multiqc && for file in "${readduplication_multiqc[@]}"; do [ -e "$file" ] && cp -r $file $multiqc_tmpdir/; done
IFS="," read -ra tin_multiqc <<< $par_tin_multiqc && for file in "${tin_multiqc[@]}"; do [ -e "$file" ] && cp -r $file $multiqc_tmpdir/; done

multiqc_config="assets/multiqc_config.yml"

multiqc \
    -f \
    ${par_multiqc_title:+--title $par_multiqc_title} \
    ${par_multiqc_custom_config:+--config $par_multiqc_custom_config} \
    $multiqc_tmpdir

# cat <<-END_VERSIONS > versions.yml
# "$key":
#     multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
# END_VERSIONS
