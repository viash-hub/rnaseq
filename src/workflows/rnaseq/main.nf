workflow run_wf {
  take:
    input_ch

  main:
    reference_ch = input_ch

    | map { id, state ->
      def biotype = state.gencode ? "gene_type" : state.featurecounts_group_type 
      def filter_gtf = 
        (
          ( !state.skip_alignment && state.aligner) // Condition 1: Alignment is required and aligner is set
          || ( !state.skip_pseudo_alignment && state.pseudo_aligner ) // Condition 2: Pseudoalignment is required and pseudoaligner is set 
          || ( !state.transcript_fasta )  // Condition 3: Transcript FASTA file is not provided
        ) 
        &&
        ( !state.skip_gtf_filter )  // Condition 4: --skip_gtf_filter is not provided
      [ id, state + [ biotype: biotype, filter_gtf: filter_gtf ] ]
    } 
    
    | toSortedList
    
    | map { list -> 
        [ "ref",  
          [ fasta: list[1][-1].fasta,
          gtf: list[1][-1].gtf, 
          gff: list[1][-1].gff, 
          additional_fasta: list[1][-1].additional_fasta,
          transcript_fasta: list[1][-1].transcript_fasta, 
          gene_bed: list[1][-1].gene_bed,
          splicesites: list[1][-1].splicesites,
          bbsplit_fasta_list: list[1][-1].bbsplit_fasta_list,
          star_index: list[1][-1].star_index,
          rsem_index: list[1][-1].rsem_index,
          salmon_index: list[1][-1].salmon_index,
          hisat2_index: list[1][-1].hisat2_index,
          bbsplit_index: list[1][-1].bbsplit_index,
          skip_bbsplit: list[1][-1].skip_bbsplit,
          gencode: list[1][-1].gencode,
          biotype: list[1][-1].biotype, 
          filter_gtf: list[1][-1].filter_gtf,
          pseudo_aligner_kmer_size: list[1][-1].pseudo_aligner_kmer_size ]
        ]
    } 

    // prepare all the necessary files for reference genome
    | prepare_genome.run ( 
        fromState: [
          "fasta": "fasta", 
          "gtf": "gtf", 
          "gff": "gff",
          "additional_fasta": "additional_fasta", 
          "transcript_fasta": "transcript_fasta", 
          "gene_bed": "gene_bed",
          "bbsplit_fasta_list": "bbsplit_fasta_list", 
          "star_index": "star_index", 
          "rsem_index": "rsem_index",
          "salmon_index": "salmon_index",
          "kallisto_index": "kallisto_index",
          "pseudo_aligner_kmer_size": "pseudo_aligner_kmer_size",
          // "splicesites": "splicesites",
          // "hisat2_index": "hisat2_index",
          "bbsplit_index": "bbsplit_index",
          "skip_bbsplit": "skip_bbsplit", 
          "gencode": "gencode", 
          "biotype": "biotype",
          "filter_gtf": "filter_gtf",
          "versions": "versions" 
        ],
        toState: [
          "fasta": "uncompressed_fasta", 
          "gtf": "gtf_uncompressed", 
          "transcript_fasta": "transcript_fasta_uncompressed", 
          "fai": "fai", 
          "chrom_sizes": "chrom_sizes", 
          "bbsplit_index": "bbsplit_index_uncompressed", 
          "star_index": "star_index_uncompressed", 
          "salmon_index": "salmon_index_uncompressed", 
          "kallisto_index": "kallisto_index_uncompressed",
          "gene_bed": "gene_bed_uncompressed",
          "versions": "updated_versions" 
        ]
    )

    // Check if contigs in genome fasta file > 512 Mbp
    | map { id, state -> 
      (isBelowMaxContigSize(state.fai)) ? [id, state] : [id, state + [bam_csi_index: true]]
    }

    | map { list -> list[1]}

    analysis_ch = input_ch

    | combine(reference_ch)

    | map { list -> [list[0], list[1] + list[2]] }

    // Concatenate FastQ files from same sample if required
    | cat_fastq.run (
        fromState: [
          "read_1": "fastq_1", 
          "read_2": "fastq_2",
          "versions": "versions"
        ], 
        toState: [ 
          "fastq_1": "fastq_1",
          "fastq_2": "fastq_2", 
          "versions": "updated_versions"
        ]
    )

    // Pre-process fastq files
    | pre_processing.run ( 
        fromState: [
          "id": "id", 
          "fastq_1": "fastq_1",
          "fastq_2": "fastq_2", 
          "umitools_bc_pattern": "umitools_bc_pattern",
          "umitools_bc_pattern2": "umitools_bc_pattern2",
          "strandedness": "strandedness",
          "transcript_fasta": "transcript_fasta", 
          "gtf": "gtf",
          "with_umi": "with_umi", 
          "bbsplit_index": "bbsplit_index", 
          "bbsplit_fasta_list": "bbsplit_fasta_list", 
          "bc_pattern": "bc_pattern", 
          "ribo_database_manifest": "ribo_database_manifest", 
          "salmon_index": "salmon_index",
          "skip_qc": "skip_qc",
          "skip_fastqc": "skip_fastqc",
          "skip_skip_umi_extract": "skip_umi_extract",
          "umi_discard_read": "umi_discard_read",
          "skip_trimming": "skip_trimming",
          "trimmer": "trimmer",
          "skip_bbsplit": "skip_bbsplit",
          "remove_ribo_rna": "remove_ribo_rna", 
          "versions": "versions"
        ], 
        toState: [ 
          "fastqc_html_1": "fastqc_html_1",
          "fastqc_html_2": "fastqc_html_2",
          "fastqc_zip_1": "fastqc_zip_1",
          "fastqc_zip_2": "fastqc_zip_2",  
          "fastq_1": "qc_output1",
          "fastq_2": "qc_output2", 
          "trim_log_1": "trim_log_1", 
          "trim_log_2": "trim_log_2", 
          "trim_zip_1": "trim_zip_1",
          "trim_zip_2": "trim_zip_2",
          "trim_html_1": "trim_html_1",
          "trim_html_2": "trim_html_2",
          "passed_trimmed_reads": "passed_trimmed_reads",
          "num_trimmed_reads": "num_trimmed_reads",
          "sortmerna_log": "sortmerna_log",
          "salmon_json_info": "salmon_json_info", 
          "fastp_failed_trim": "failed_trim",
          "fastp_failed_trim_unpaired1": "failed_trim_unpaired1",
          "fastp_failed_trim_unpaired2": "failed_trim_unpaired2",
          "fastp_trim_json": "trim_json",
          "fastp_trim_html": "trim_html",
          "fastp_trim_merged_out": "trim_merged_out",
          "versions": "updated_versions"
        ]
    )

    // Infer strandedness from Salmon pseudo-alignment results
    | map { id, state -> 
    (state.strandedness == 'auto') ? 
      [ id, state + [strandedness: getSalmonInferredStrandedness(state.salmon_json_info)] ] : 
      [id, state] 
    }

    // Filter FastQ files based on minimum trimmed read count after adapter trimming
    | map { id, state -> 
      def input = state.fastq_2 ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
      def num_reads = (state.skip_trimming) ? 
        state.min_trimmed_reads + 1 : 
        (
          (state.trimmer == "fastp") ? 
            getFastpReadsAfterFiltering(state.fastp_trim_json) : 
            (
              (!state.skip_trimming && input.size() == 2) ?
                getTrimGaloreReadsAfterFiltering(state.trim_log_2) : 
                getTrimGaloreReadsAfterFiltering(state.trim_log_1)
            )
        )
      def passed_trimmed_reads = 
        (state.skip_trimming || (num_reads >= state.min_trimmed_reads)) ? 
          true : 
          false 
      [ id, state + [num_trimmed_reads: num_reads, passed_trimmed_reads: passed_trimmed_reads] ] 
    }

    // Genome alignment and quantification
    | genome_alignment_and_quant.run (
        runIf: { id, state -> !state.skip_alignment && state.passed_trimmed_reads },
        fromState: [
          "id": "id", 
          "fastq_1": "fastq_1",
          "fastq_2": "fastq_2", 
          "strandedness": "strandedness", 
          "gtf": "gtf",
          "transcript_fasta": "transcript_fasta",
          "bam_csi_index": "bam_csi_index", 
          "aligner": "aligner",
          "rsem_index": "rsem_index",
          "star_index": "star_index", 
          "extra_star_align_args": "extra_star_align_args", 
          "star_ignore_sjdbgtf": "star_ignore_sjdbgtf",
          "seq_platform": "seq_platform", 
          "seq_center": "seq_center",
          "with_umi": "with_umi", 
          "umi_dedup_stats": "umi_dedup_stats",
          "gtf_group_features": "gtf_group_features",
          "gtf_extra_attributes": "gtf_extra_attributes",
          "salmon_quant_libtype": "salmon_quant_libtype", 
          "versions": "versions" 
        ],
        toState: [
          "star_alignment": "star_alignment", 
          "star_multiqc": "star_multiqc", 
          "genome_bam_sorted": "genome_bam_sorted",
          "genome_bam_index": "genome_bam_index", 
          "genome_bam_stats": "genome_bam_stats", 
          "genome_bam_flagstat": "genome_bam_flagstat", 
          "genome_bam_idxstats": "genome_bam_idxstats", 
          "transcriptome_bam": "transcriptome_bam", 
          "transcriptome_bam_index": "transcriptome_bam_index", 
          "transcriptome_bam_stats": "transcriptome_bam_stats", 
          "transcriptome_bam_flagstat": "transcriptome_bam_flagstat", 
          "transcriptome_bam_idxstats": "transcriptome_bam_idxstats",
          "quant_results": "quant_results", 
          "versions": "updated_versions"
        ]
    )

    // Filter channels to get samples that passed STAR minimum mapping percentage
    | map { id, state -> 
      def percent_mapped = getStarPercentMapped(state.star_multiqc) 
      def passed_mapping = (percent_mapped >= state.min_mapped_reads) ? true : false
      [ id, state + [percent_mapped: percent_mapped, passed_mapping: passed_mapping] ]
    }

    // Pseudo-alignment and quantification
    | pseudo_alignment_and_quant.run (
        runIf: { id, state -> !state.skip_pseudo_alignment && state.passed_trimmed_reads },
        fromState: [
          "id": "id", 
          "fastq_1": "fastq_1",
          "fastq_2": "fastq_2", 
          "strandedness": "strandedness", 
          "gtf": "gtf",
          "transcript_fasta": "transcript_fasta",
          "pseudo_aligner": "pseudo_aligner",
          "pseudo_index": "pseudo_index",
          "extra_star_align_args": "extra_star_align_args", 
          "star_ignore_sjdbgtf": "star_ignore_sjdbgtf",
          "seq_platform": "seq_platform", 
          "seq_center": "seq_center",
          "with_umi": "with_umi", 
          "umi_dedup_stats": "umi_dedup_stats",
          "gtf_group_features": "gtf_group_features",
          "gtf_extra_attributes": "gtf_extra_attributes",
          "libtype": "salmon_quant_libtype", 
          "kallisto_quant_fragment_length": "kallisto_quant_fragment_length",
          "kallisto_quant_fragment_length_sd": "kallisto_quant_fragment_length_sd",
          "versions": "versions" 
        ],
        toState: [
          "pseudo_multiqc": "pseudo_multiqc", 
          "quant_results": "quant_results", 
          "versions": "updated_versions"
        ]
    )

    | map { id, state ->
      def input = state.fastq_2 ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
      def paired = input.size() == 2
      [ id, state + [ paired: paired ] ]
    }

    // Post-processing
    | post_processing.run (
        runIf: { id, state -> !state.skip_alignment && state.passed_trimmed_reads && state.passed_mapping },
        fromState: [
          "id": "id", 
          "paired": "paired", 
          "strandedness": "strandedness", 
          "fasta": "fasta",
          "fai": "fai", 
          "gtf": "gtf", 
          "genome_bam": "genome_bam_sorted", 
          "chrom_sizes": "chrom_sizes", 
          "star_multiqc": "star_multiqc",
          "extra_picard_args": "extra_picard_args", 
          "extra_stringtie_args": "extra_stringtie_args", 
          "stringtie_ignore_gtf": "stringtie_ignore_gtf", 
          "extra_bedtools_args": "extra_bedtools_args", 
          "bam_csi_index": "bam_csi_index", 
          "min_mapped_reads": "min_mapped_reads", 
          "with_umi": "with_umi",
          "skip_qc": "skip_qc",
          "skip_markduplicates": "skip_markduplicates", 
          "skip_stringtie": "skip_stringtie", 
          "skip_bigwig":"gencode", 
          "versions": "versions"
        ], 
        toState: [
          "genome_bam_sorted": "processed_genome_bam", 
          "genome_bam_index": "genome_bam_index",
          "genome_bam_stats": "genome_bam_stats",
          "genome_bam_flagstat": "genome_bam_flagstat", 
          "genome_bam_idxstats": "genome_bam_idxstats", 
          "markduplicates_metrics": "markduplicates_metrics",
          "stringtie_transcript_gtf": "stringtie_transcript_gtf",
          "stringtie_coverage_gtf": "stringtie_coverage_gtf",
          "stringtie_abundance": "stringtie_abundance",
          "stringtie_ballgown": "stringtie_ballgown", 
          "bedgraph_forward": "bedgraph_forward",
          "bedgraph_reverse": "bedgraph_reverse",
          "bigwig_forward": "bigwig_forward",
          "bigwig_reverse": "bigwig_reverse", 
          "versions": "updated_versions"
        ]
    )

    // Final QC
    | quality_control.run (
        fromState: [
          "id": "id", 
          "paired": "paired", 
          "strandedness": "strandedness", 
          "skip_align": "skip_alignment",
          "skip_pseudo_align": "skippseudo_alignment",
          "gtf": "gtf", 
          "genome_bam": "genome_bam_sorted", 
          "genome_bam_index": "genome_bam_index",
          "quant_results": "quant_results", 
          "gene_bed": "gene_bed",
          "extra_preseq_args": "extra_preseq_args",
          "extra_featurecounts_args": "extra_featurecounts_args", 
          "biotype": "biotype", 
          "biotypes_header": "biotypes_header",
          "skip_biotype_qc": "skip_biotype_qc", 
          "featurecounts_group_type": "featurecounts_group_type", 
          "featurecounts_feature_type": "featurecounts_feature_type", 
          "gencode": "gencode",
          "skip_deseq2_qc": "skip_deseq2_qc",  
          "extra_deseq2_args": "extra_deseq2_args",
          "extra_deseq2_args2": "extra_deseq2_args2",
          "multiqc_custom_config": "multiqc_custom_config", 
          "multiqc_title": "multiqc_title", 
          "multiqc_methods_description": "multiqc_methods_description",
          "fastqc_zip_1": "fastqc_zip_1",
          "fastqc_zip_2": "fastqc_zip_2",  
          "trim_log_1": "trim_log_1", 
          "trim_log_2": "trim_log_2", 
          "trim_zip_1": "trim_zip_1",
          "trim_zip_2": "trim_zip_2",
          "sortmerna_multiqc": "sortmerna_log", 
          "star_multiqc": "star_multiqc", 
          "genome_bam_stats": "genome_bam_stats", 
          "genome_bam_flagstat": "genome_bam_flagstat", 
          "genome_bam_idxstats": "genome_bam_idxstats", 
          "markduplicates_multiqc": "markduplicates_metrics", 
          "rseqc_modules": "rseqc_modules",
          "num_trimmed_reads": "num_trimmed_reads",
          "passed_trimmed_reads": "passed_trimmed_reads",
          "passed_mapping": "passed_mapping",
          "percent_mapped": "percent_mapped", 
          "versions": "versions"
        ], 
        toState: [
          "preseq_output": "preseq_output",
          "bamstat_output": "bamstat_output",
          "strandedness_output": "strandedness_output",
          "inner_dist_output_stats": "inner_dist_output_stats",
          "inner_dist_output_dist": "inner_dist_output_dist",
          "inner_dist_output_freq": "inner_dist_output_freq",
          "inner_dist_output_plot": "inner_dist_output_plot",
          "inner_dist_output_plot_r": "inner_dist_output_plot_r",
          "junction_annotation_output_log": "junction_annotation_output_log",
          "junction_annotation_output_plot_r": "junction_annotation_output_plot_r",
          "junction_annotation_output_junction_bed": "junction_annotation_output_junction_bed",
          "junction_annotation_output_junction_interact": "junction_annotation_output_junction_interact",
          "junction_annotation_output_junction_sheet": "junction_annotation_output_junction_sheet",
          "junction_annotation_output_splice_events_plot": "junction_annotation_output_splice_events_plot",
          "junction_annotation_output_splice_junctions_plot": "junction_annotation_output_splice_junctions_plot",
          "junction_saturation_output_plot_r": "junction_saturation_output_plot_r",
          "junction_saturation_output_plot": "junction_saturation_output_plot",
          "read_distribution_output": "read_distribution_output",
          "read_duplication_output_duplication_rate_plot_r": "read_duplication_output_duplication_rate_plot_r",
          "read_duplication_output_duplication_rate_plot": "read_duplication_output_duplication_rate_plot",
          "read_duplication_output_duplication_rate_mapping": "read_duplication_output_duplication_rate_mapping",
          "read_duplication_output_duplication_rate_sequence": "read_duplication_output_duplication_rate_sequence",
          "tin_output_summary": "tin_output_summary",
          "tin_output_metrics": "tin_output_metrics",
          "dupradar_output_dupmatrix": "dupradar_output_dupmatrix",
          "dupradar_output_dup_intercept_mqc": "dupradar_output_dup_intercept_mqc",
          "dupradar_output_duprate_exp_boxplot": "dupradar_output_duprate_exp_boxplot",
          "dupradar_output_duprate_exp_densplot": "dupradar_output_duprate_exp_densplot",
          "dupradar_output_duprate_exp_denscurve_mqc": "dupradar_output_duprate_exp_denscurve_mqc",
          "dupradar_output_expression_histogram": "dupradar_output_expression_histogram",
          "dupradar_output_intercept_slope": "dupradar_output_intercept_slope",
          "qualimap_output_dir": "qualimap_output_dir",
          "qualimap_output_pdf": "qualimap_output_pdf",
          "featurecounts": "featurecounts",
          "featurecounts_summary": "featurecounts_summary",
          "featurecounts_multiqc": "featurecounts_multiqc", 
          "featurecounts_rrna_multiqc": "featurecounts_rrna_multiqc",
          "tpm_gene": "tpm_gene",
          "counts_gene": "counts_gene",
          "counts_gene_length_scaled": "counts_gene_length_scaled",
          "counts_gene_scaled": "counts_gene_scaled", 
          "tpm_transcript": "tpm_transcript", 
          "counts_transcript": "counts_transcript", 
          "salmon_merged_summarizedexperiment": "salmon_merged_summarizedexperiment",
          "deseq2_output": "deseq2_output", 
          "multiqc_report": "multiqc_report", 
          "multiqc_data": "multiqc_data", 
          "multiqc_plots": "multiqc_plots",
          "multiqc_versions": "multiqc_versions", 
          "versions": "updated_versions" 
        ] 
    )

    | map { id, state -> 
      def mod_state = state.findAll { key, value -> value instanceof java.nio.file.Path && value.exists() }
      [ id, mod_state ]
    }

    | setState (
      [
        "output_fasta": "fasta", 
        "output_gtf": "gtf", 
        "output_transcript_fasta": "transcript_fasta", 
        "output_gene_bed": "gene_bed", 
        "output_bbsplit_index": "bbsplit_index", 
        "output_star_index": "star_index", 
        "output_salmon_index": "salmon_index",
        "output_kallisto_index": "kallisto_index",
        "fastqc_html_1": "fastqc_html_1",
        "fastqc_html_2": "fastqc_html_2",
        "fastqc_zip_1": "fastqc_zip_1",
        "fastqc_zip_2": "fastqc_zip_2",  
        "output_fastq_1": "fastq_1",
        "output_fastq_2": "fastq_2", 
        "trim_log_1": "trim_log_1", 
        "trim_log_2": "trim_log_2", 
        "trim_zip_1": "trim_zip_1",
        "trim_zip_2": "trim_zip_2",
        "trim_html_1": "trim_html_1",
        "trim_html_2": "trim_html_2",
        "sortmerna_log": "sortmerna_log",
        "star_alignment": "star_alignment", 
        "genome_bam_sorted": "genome_bam_sorted",
        "genome_bam_index": "genome_bam_index", 
        "genome_bam_stats": "samtools_stats", 
        "genome_bam_flagstat": "samtools_flagstat", 
        "genome_bam_idxstats": "samtools_idxstats", 
        "transcriptome_bam": "transcriptome_bam", 
        "transcriptome_bam_index": "transcriptome_bam_index", 
        "transcriptome_bam_stats": "transcriptome_bam_stats", 
        "transcriptome_bam_flagstat": "transcriptome_bam_flagstat", 
        "transcriptome_bam_idxstats": "transcriptome_bam_idxstats",
        "salmon_quant_results": "salmon_quant_results",
        "stringtie_transcript_gtf": "stringtie_transcript_gtf",
        "stringtie_coverage_gtf": "stringtie_coverage_gtf",
        "stringtie_abundance": "stringtie_abundance",
        "stringtie_ballgown": "stringtie_ballgown", 
        "featurecounts": "featurecounts",
        "featurecounts_summary": "featurecounts_summary", 
        "featurecounts_multiqc": "featurecounts_multiqc", 
        "featurecounts_rrna_multiqc": "featurecounts_rrna_multiqc", 
        "bedgraph_forward": "bedgraph_forward",
        "bedgraph_reverse": "bedgraph_reverse",
        "bigwig_forward": "bigwig_forward",
        "bigwig_reverse": "bigwig_reverse",
        "preseq_output": "preseq_output",
        "bamstat_output": "bamstat_output",
        "strandedness_output": "strandedness_output",
        "inner_dist_output_stats": "inner_dist_output_stats",
        "inner_dist_output_dist": "inner_dist_output_dist",
        "inner_dist_output_freq": "inner_dist_output_freq",
        "inner_dist_output_plot": "inner_dist_output_plot",
        "inner_dist_output_plot_r": "inner_dist_output_plot_r",
        "junction_annotation_output_log": "junction_annotation_output_log",
        "junction_annotation_output_plot_r": "junction_annotation_output_plot_r",
        "junction_annotation_output_junction_bed": "junction_annotation_output_junction_bed",
        "junction_annotation_output_junction_interact": "junction_annotation_output_junction_interact",
        "junction_annotation_output_junction_sheet": "junction_annotation_output_junction_sheet",
        "junction_annotation_output_splice_events_plot": "junction_annotation_output_splice_events_plot",
        "junction_annotation_output_splice_junctions_plot": "junction_annotation_output_splice_junctions_plot",
        "junction_saturation_output_plot_r": "junction_saturation_output_plot_r",
        "junction_saturation_output_plot": "junction_saturation_output_plot",
        "read_distribution_output": "read_distribution_output",
        "read_duplication_output_duplication_rate_plot_r": "read_duplication_output_duplication_rate_plot_r",
        "read_duplication_output_duplication_rate_plot": "read_duplication_output_duplication_rate_plot",
        "read_duplication_output_duplication_rate_mapping": "read_duplication_output_duplication_rate_mapping",
        "read_duplication_output_duplication_rate_sequence": "read_duplication_output_duplication_rate_sequence",
        "tin_output_summary": "tin_output_summary",
        "tin_output_metrics": "tin_output_metrics",
        "dupradar_output_dupmatrix": "dupradar_output_dupmatrix",
        "dupradar_output_dup_intercept_mqc": "dupradar_output_dup_intercept_mqc",
        "dupradar_output_duprate_exp_boxplot": "dupradar_output_duprate_exp_boxplot",
        "dupradar_output_duprate_exp_densplot": "dupradar_output_duprate_exp_densplot",
        "dupradar_output_duprate_exp_denscurve_mqc": "dupradar_output_duprate_exp_denscurve_mqc",
        "dupradar_output_expression_histogram": "dupradar_output_expression_histogram",
        "dupradar_output_intercept_slope": "dupradar_output_intercept_slope",
        "qualimap_output_dir": "qualimap_output_dir",
        "qualimap_output_pdf": "qualimap_output_pdf", 
        "tpm_gene": "tpm_gene",
        "counts_gene": "counts_gene",
        "counts_gene_length_scaled": "counts_gene_length_scaled",
        "counts_gene_scaled": "counts_gene_scaled", 
        "tpm_transcript": "tpm_transcript", 
        "counts_transcript": "counts_transcript", 
        "salmon_merged_summarizedexperiment": "salmon_merged_summarizedexperiment",
        "deseq2_output": "deseq2_output", 
        "multiqc_report": "multiqc_report", 
        "multiqc_data": "multiqc_data", 
        "multiqc_plots": "multiqc_plots",
        "multiqc_versions": "multiqc_versions", 
        "versions": "versions"
      ]
    )

    output_ch = analysis_ch

  emit:
    output_ch
}

import nextflow.Nextflow
//
// Function to generate an error if contigs in genome fasta file > 512 Mbp
//
def isBelowMaxContigSize(fai_file) {
  def max_size = 512000000
  fai_file.eachLine { line ->
    def lspl  = line.split('\t')
    def chrom = lspl[0]
    def size  = lspl[1]
    if (size.toInteger() > max_size) {
      def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
        "  Contig longer than ${max_size}bp found in reference genome!\n\n" +
        "  ${chrom}: ${size}\n\n" +
        "  Provide the '--bam_csi_index' parameter to use a CSI instead of BAI index.\n\n"
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
      Nextflow.error(error_string)
      return false
    }
  }
  return true
}

import groovy.json.JsonSlurper
//
// Function that parses Salmon quant 'meta_info.json' output file to get inferred strandedness
//
def getSalmonInferredStrandedness(json_file) {
  def lib_type = new JsonSlurper().parseText(json_file.text).get('library_types')[0]
  def strandedness = 'reverse'
  if (lib_type) {
    if (lib_type in ['U', 'IU']) {
      strandedness = 'unstranded'
    } 
    else if (lib_type in ['SF', 'ISF']) {
      strandedness = 'forward'
    } 
    else if (lib_type in ['SR', 'ISR']) {
      strandedness = 'reverse'
    }
  }
  return strandedness
}

//
// Function that parses TrimGalore log output file to get total number of reads after trimming
//
def getTrimGaloreReadsAfterFiltering(log_file) {
  def total_reads = 0
  def filtered_reads = 0
  log_file.eachLine { line ->
    def total_reads_matcher = line =~ /([\d\.]+)\ssequences processed in total/
    def filtered_reads_matcher = line =~ /shorter than the length cutoff[^:]+:\s([\d\.]+)/
    if (total_reads_matcher) total_reads = total_reads_matcher[0][1].toFloat()
    if (filtered_reads_matcher) filtered_reads = filtered_reads_matcher[0][1].toFloat()
  }
  return total_reads - filtered_reads
}

//
// Function that parses fastp json output file to get total number of reads after trimming
//
def getFastpReadsAfterFiltering(json_file) {
    def Map json = (Map) new JsonSlurper().parseText(json_file.text).get('summary')
    return json['after_filtering']['total_reads'].toLong()
}

//
// Function that parses and returns the alignment rate from the STAR log output
//
def getStarPercentMapped(align_log) {
  def percent_aligned = 0
  def pattern = /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/
  align_log.eachLine { line ->
    def matcher = line =~ pattern
    if (matcher) {
        percent_aligned = matcher[0][1].toFloat()
    }
  }
  return percent_aligned
}
