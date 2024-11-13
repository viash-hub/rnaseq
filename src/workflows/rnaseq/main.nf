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
          [ fasta: list.collect { id, state -> state.fasta }.unique()[0],
          gtf: list.collect { id, state -> state.gtf }.unique()[0], 
          gff: list.collect { id, state -> state.gff }.unique()[0], 
          additional_fasta: list.collect { id, state -> state.additional_fasta }.unique()[0],
          transcript_fasta:list.collect { id, state -> state.transcript_fasta }.unique()[0], 
          gene_bed: list.collect { id, state -> state.gene_bed }.unique()[0],
          bbsplit_fasta_list: list.collect { id, state -> state.bbsplit_fasta_list }.unique()[0],
          aligner: list.collect { id, state -> state.aligner }.unique()[0],
          pseudo_aligner: list.collect { id, state -> state.pseudo_aligner }.unique()[0],
          star_index: list.collect { id, state -> state.star_index }.unique()[0],
          rsem_index: list.collect { id, state -> state.rsem_index }.unique()[0],
          salmon_index: list.collect { id, state -> state.salmon_index }.unique()[0],
          kallisto_index: list.collect { id, state -> state.kallisto_index }.unique()[0],
          // splicesites: list.collect { id, state -> state.splicesites }.unique()[0],
          // hisat2_index: list.collect { id, state -> state.hisat2_index }.unique()[0],
          bbsplit_index: list.collect { id, state -> state.bbsplit_index }.unique()[0],
          skip_bbsplit: list.collect { id, state -> state.skip_bbsplit }.unique()[0],
          gencode: list.collect { id, state -> state.gencode }.unique()[0],
          biotype: list.collect { id, state -> state.biotype }.unique()[0], 
          filter_gtf: list.collect { id, state -> state.filter_gtf }.unique()[0],
          pseudo_aligner_kmer_size: list.collect { id, state -> state.pseudo_aligner_kmer_size }.unique()[0] ]
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
          "aligner": "aligner",
          "pseudo_aligner": "pseudo_aligner",
          "skip_alignment": "skip_alignment"
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
          "gene_bed": "gene_bed_uncompressed"
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
          "read_2": "fastq_2"
        ], 
        toState: [ 
          "fastq_1": "fastq_1",
          "fastq_2": "fastq_2"
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
          "remove_ribo_rna": "remove_ribo_rna"
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
          "salmon_quant_output": "salmon_quant_output",
          "fastp_failed_trim": "failed_trim",
          "fastp_failed_trim_unpaired1": "failed_trim_unpaired1",
          "fastp_failed_trim_unpaired2": "failed_trim_unpaired2",
          "fastp_trim_json": "trim_json",
          "fastp_trim_html": "trim_html",
          "fastp_trim_merged_out": "trim_merged_out"
        ]
    )

    // Infer strandedness from Salmon pseudo-alignment results
    | map { id, state -> 
    (state.strandedness == 'auto') ? 
      [ id, state + [strandedness: getSalmonInferredStrandedness(state.salmon_quant_output)] ] : 
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
          "salmon_index": "salmon_index",
          "extra_rsem_calculate_expression_args": "extra_rsem_calculate_expression_args" 
        ],
        toState: [
          "star_multiqc": "star_multiqc", 
          "rsem_multiqc": "rsem_multiqc",
          "salmon_multiqc": "salmon_multiqc",
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
          "quant_out_dir": "quant_out_dir",
          "quant_results_file": "quant_results_file",
          "rsem_counts_gene": "rsem_counts_gene",
          "rsem_counts_transcripts": "rsem_counts_transcripts",
          "bam_genome_rsem": "bam_genome_rsem",
          "bam_transcript_rsem": "bam_transcript_rsem"
        ]
    )

    // Filter channels to get samples that passed STAR minimum mapping percentage
    | map { id, state -> 
      def percent_mapped = (!state.skip_alignment) ? getStarPercentMapped(state.star_multiqc) : 0.0
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
          "salmon_index": "salmon_index",
          "kallisto_index": "kallisto_index",
          "extra_star_align_args": "extra_star_align_args", 
          "star_ignore_sjdbgtf": "star_ignore_sjdbgtf",
          "seq_platform": "seq_platform", 
          "seq_center": "seq_center",
          "with_umi": "with_umi", 
          "umi_dedup_stats": "umi_dedup_stats",
          "gtf_group_features": "gtf_group_features",
          "gtf_extra_attributes": "gtf_extra_attributes",
          "lib_type": "salmon_quant_libtype", 
          "kallisto_quant_fragment_length": "kallisto_quant_fragment_length",
          "kallisto_quant_fragment_length_sd": "kallisto_quant_fragment_length_sd"
        ],
        toState: [
          "pseudo_quant_out_dir": "quant_out_dir",
          "pseudo_salmon_quant_results_file": "salmon_quant_results_file",
          "pseudo_kallisto_quant_results_file": "kallisto_quant_results_file",
          "pseudo_multiqc": "pseudo_multiqc"
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
          "skip_bigwig":"gencode"
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
          "bigwig_reverse": "bigwig_reverse"
        ]
    )

    // Final QC
    | quality_control.run (
        fromState: [
          "id": "id", 
          "paired": "paired", 
          "strandedness": "strandedness", 
          "skip_align": "skip_alignment",
          "skip_pseudo_align": "skip_pseudo_alignment",
          "skip_dupradar": "skip_dupradar",
          "skip_qualimap": "skip_qualimap",
          "skip_rseqc": "skip_rseqc",
          "skip_multiqc": "skip_multiqc",
          "skip_preseq": "skip_preseq",
          "gtf": "gtf", 
          "num_trimmed_reads": "num_trimmed_reads",
          "passed_trimmed_reads": "passed_trimmed_reads",
          "passed_mapping": "passed_mapping",
          "percent_mapped": "percent_mapped",
          "genome_bam": "genome_bam_sorted", 
          "genome_bam_index": "genome_bam_index",
          "salmon_multiqc": "salmon_multiqc",
          "quant_results_file": "quant_results_file",
          "rsem_multiqc": "rsem_multiqc",
          "rsem_counts_gene": "rsem_counts_gene",
          "rsem_counts_transcripts": "rsem_counts_transcripts",
          "pseudo_multiqc": "pseudo_multiqc",
          "pseudo_quant_out_dir": "pseudo_quant_out_dir",
          "pseudo_salmon_quant_results_file": "pseudo_salmon_quant_results_file", 
          "pseudo_kallisto_quant_results_file": "pseudo_kallisto_quant_results_file", 
          "aligner": "aligner",
          "pseudo_aligner": "pseudo_aligner",
          "gene_bed": "gene_bed",
          "extra_preseq_args": "extra_preseq_args",
          "biotype": "biotype", 
          "skip_biotype_qc": "skip_biotype_qc", 
          "featurecounts_group_type": "featurecounts_group_type", 
          "featurecounts_feature_type": "featurecounts_feature_type", 
          "gencode": "gencode",
          "skip_deseq2_qc": "skip_deseq2_qc",  
          "deseq2_vst": "deseq2_vst",
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
          "rseqc_modules": "rseqc_modules"
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
          "qunat_merged_summarizedexperiment": "quant_merged_summarizedexperiment",
          "deseq2_output": "deseq2_output", 
          "multiqc_report": "multiqc_report", 
          "multiqc_data": "multiqc_data", 
          "multiqc_plots": "multiqc_plots"
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
        "star_log": "star_multiqc", 
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
        "salmon_quant_results": "quant_out_dir",
        "pseudo_quant_results": "pseudo_quant_out_dir",
        "markduplicates_metrics": "markduplicates_metrics",
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
        "quant_merged_summarizedexperiment": "quant_merged_summarizedexperiment",
        "deseq2_output": "deseq2_output", 
        "pseudo_tpm_gene": "pseudo_tpm_gene",
        "pseudo_counts_gene": "pseudo_counts_gene",
        "pseudo_counts_gene_length_scaled": "pseudo_counts_gene_length_scaled",
        "pseudo_counts_gene_scaled": "pseudo_counts_gene_scaled", 
        "pseudo_tpm_transcript": "pseudo_tpm_transcript", 
        "pseudo_counts_transcript": "pseudo_counts_transcript", 
        "pseudo_lengths_gene": "pseudo_lengths_gene",
        "pseudo_lengths_transcript": "pseudo_lengths_transcript",
        "pseudo_quant_merged_summarizedexperiment": "pseudo_quant_merged_summarizedexperiment",
        "deseq2_output_pseudo": "deseq2_output_pseudo",  
        "multiqc_report": "multiqc_report", 
        "multiqc_data": "multiqc_data", 
        "multiqc_plots": "multiqc_plots"
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
def getSalmonInferredStrandedness(salmon_quant_output) {
  def json_file = new File(salmon_quant_output).listFiles().find { it.name == "meta_info.json" || it.isDirectory() && it.listFiles().find { it.name == "meta_info.json" } }
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
// Function that parses and returns the alignment rate from the STAR log outputs
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
