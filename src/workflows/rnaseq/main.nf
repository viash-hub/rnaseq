workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch

    | map { id, state ->
      def biotype = state.gencode ? "gene_type" : state.featurecounts_group_type 
      def bc_pattern = state.umitools_bc_pattern2 ? [ state.umitools_bc_pattern, state.umitools_bc_pattern2 ] : [ state.umitools_bc_pattern ] 
      [ id, state + [ biotype: biotype, bc_pattern: bc_pattern ] ]
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
          "splicesites": "splicesites",
          "bbsplit_fasta_list": "bbsplit_fasta_list", 
          "star_index": "star_index", 
          "rsem_index": "rsem_index",
          "salmon_index": "salmon_index",
          "hisat2_index": "hisat2_index",
          "bbsplit_index": "bbsplit_index",
          "gencode": "gencode", 
          "biotype": "biotype" 
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
          "gene_bed": "gene_bed_uncompressed" 
        ]
    )

    // Check if contigs in genome fasta file > 512 Mbp
    | map { id, state -> 
      (isBelowMaxContigSize(state.fai)) ? [id, state] : [id, state + [bam_csi_index: true]]
    }

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
          "salmon_index": "salmon_index"
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
          "salmon_json_info": "salmon_json_info"
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
      def num_reads = state.min_trimmed_reads + 1
      num_reads = 
        (!state.skip_trimming && input.size() == 2) ?
          getTrimGaloreReadsAfterFiltering(state.trim_log_2) : 
          getTrimGaloreReadsAfterFiltering(state.trim_log_1)
      def passed_trimmed_reads = 
        (state.skip_trimming || (num_reads >= state.min_trimmed_reads)) ? 
          true : 
          false 
      [ id, state + [num_trimmed_reads: num_reads, passed_trimmed_reads: passed_trimmed_reads] ] 
    }
    // | filter { id, state -> state.skip_trimming || state.passed_trimmed_reads }
    // TODO: Get list of samples that failed trimming threshold for MultiQC report

    // Genome alignment and quantification
    | genome_alignment_and_quant.run (
        fromState: [
          "id": "id", 
          "fastq_1": "fastq_1",
          "fastq_2": "fastq_2", 
          "strandedness": "strandedness", 
          "gtf": "gtf",
          "transcript_fasta": "transcript_fasta",
          "bam_csi_index": "bam_csi_index", 
          "star_index": "star_index", 
          "extra_star_align_args": "extra_star_align_args", 
          "star_ignore_sjdbgtf": "star_ignore_sjdbgtf",
          "seq_platform": "seq_platform", 
          "seq_center": "seq_center",
          "with_umi": "with_umi", 
          "umi_dedup_stats": "umi_dedup_stats",
          "gtf_group_features": "gtf_group_features",
          "gtf_extra_attributes": "gtf_extra_attributes",
          "salmon_quant_libtype": "salmon_quant_libtype" 
        ],
        toState: [
          "star_alignment": "star_alignment", 
          "star_multiqc": "star_multiqc", 
          "genome_bam_sorted": "genome_bam_sorted",
          "genome_bam_index": "genome_bam_index", 
          "genome_bam_stats": "genome_bam_stats", 
          "genome_bam_flagstat": "genome_bam_flagstat", 
          "genome_bam_idxstats": "genome_bam_idxstats", 
          "transcriptome_bam_sorted": "transcriptome_bam_sorted", 
          "transcriptome_bam_index": "transcriptome_bam_index", 
          "transcriptome_bam_stats": "transcriptome_bam_stats", 
          "transcriptome_bam_flagstat": "transcriptome_bam_flagstat", 
          "transcriptome_bam_idxstats": "transcriptome_bam_idxstats",
          "salmon_quant_results": "salmon_quant_results"
        ]
    )

    // Filter channels to get samples that passed STAR minimum mapping percentage
    | map { id, state -> 
      def percent_mapped = getStarPercentMapped(state.star_multiqc) 
      def passed_mapping = (percent_mapped >= state.min_mapped_reads) ? true : false
      [ id, state + [percent_mapped: percent_mapped, passed_mapping: passed_mapping] ]
    }
    // | filter { id, state -> state.passed_mapping) }
    // TODO: Get list of samples that failed mapping for MultiQC report

    | map { id, state ->
      def input = state.fastq_2 ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
      def paired = input.size() == 2
      [ id, state + [ paired: paired ] ]
    }
    // Post-processing
    | post_processing.run (
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
          "extra_featurecounts_args": "extra_featurecounts_args", 
          "bam_csi_index": "bam_csi_index", 
          "min_mapped_reads": "min_mapped_reads", 
          "with_umi": "with_umi",
          "biotype": "biotype", 
          "biotypes_header": "biotypes_header",
          "skip_qc": "skip_qc",
          "skip_markdupkicates": "skip_markdupkicates", 
          "skip_stringtie": "skip_stringtie", 
          "skip_biotype_qc": "skip_biotype_qc", 
          "skip_bigwig": "skip_bigwig", 
          "featurecounts_group_type": "featurecounts_group_type", 
          "featurecounts_feature_type": "featurecounts_feature_type", 
          "gencode": "gencode"
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
          "featurecounts": "featurecounts",
          "featurecounts_summary": "featurecounts_summary", 
          "featurecounts_multiqc": "featurecounts_multiqc", 
          "bedgraph_forward": "bedgraph_forward",
          "bedgraph_reverse": "bedgraph_reverse",
          "bigwig_forward": "bigwig_forward",
          "bigwig_reverse": "bigwig_reverse"
        ], 
    )

    // Final QC
    | quality_control.run (
        fromState: [
          "id": "id", 
          "paired": "paired", 
          "strandedness": "strandedness", 
          "gtf": "gtf", 
          "genome_bam": "genome_bam_sorted", 
          "genome_bam_index": "genome_bam_index",
          "salmon_quant_results": "salmon_quant_results", 
          "gene_bed": "gene_bed",
          "extra_preseq_args": "extra_preseq_args",
          // "skip_deseq2_qc": "skip_deseq2_qc",  
          // "extra_deseq2_args": "extra_deseq2_args",
          // "extra_deseq2_args2": "extra_deseq2_args2",
          // "multiqc_custom_config": "multiqc_custom_config", 
          // "multiqc_title": "multiqc_title", 
          // "multiqc_logo": "multiqc_logo",
          // "multiqc_methods_description": "multiqc_methods_description",
          // "fail_trimming_multiqc": "fail_trimming_multiqc", 
          // "fail_mapping_multiqc": "fail_mapping_multiqc", 
          // "fastqc_zip_1": "fastqc_zip_1",
          // "fastqc_zip_2": "fastqc_zip_2",  
          // "trim_log_1": "trim_log_1", 
          // "trim_log_2": "trim_log_2", 
          // "trim_zip_1": "trim_zip_1",
          // "trim_zip_2": "trim_zip_2",
          // "sortmerna_multiqc": "sortmerna_log", 
          // "star_multiqc": "star_multiqc", 
          // "genome_bam_stats": "genome_bam_stats", 
          // "genome_bam_flagstat": "genome_bam_flagstat", 
          // "genome_bam_idxstats": "genome_bam_idxstats", 
          // "markduplicates_multiqc": "markduplicates_metrics", 
          // "featurecounts_multiqc": "featurecounts_multiqc"        
        ], 
        toState: [
          // "multiqc_report": "multiqc_report", 
          // "multiqc_data": "multiqc_data",
          // "multiqc_plots": "multiqc_plots",
          // "multiqc_versions": "multiqc_versions",
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
          "qualimap_output_pdf": "qualimap_output_pdf"
        ] 
    )

    // | map { id, state -> 
    //   def paired_state = (!state.paired) ? 
    //     [trim_log_2: state.remove(state.trim_log_2), trim_zip_2: state.remove(state.trim_zip_2), trim_html_2: state.remove(state.trim_html_2)] : 
    //     []
    //   def qc_state = (state.skip_qc || state.skip_fastqc) ? 
    //     [fastqc_html_1: state.remove(state.fastqc_html_1), fastqc_html_2: state.remove(state.fastqc_html_2), fastqc_zip_1: state.remove(state.fastqc_zip_1), fastqc_zip_2: state.remove(state.fastqc_zip_2)] : 
    //     []
    //   def trimming_state = (state.skip_trimming) ? 
    //     [trim_html_1: state.remove(state.trim_html_1), trim_html_2: state.remove(state.trim_html_2), trim_zip_1: state.remove(state.trim_zip_1), trim_zip_2: state.remove(state.trim_zip_2), trim_log_1: state.remove(state.trim_log_1), trim_log_2: state.remove(state.trim_log_2)] : 
    //     []
    //   def sortmerna_state = (!state.remove_ribo_rna) ? [sortmerna_log: state.remove(state.sortmerna_log)] : []
    //   [ id, state + paired_state + qc_state + trimming_state + sortmerna_state ]
    // }
    // qc_ch = analysis_ch
    | toSortedList
    | map { list -> 
        def fastqc_html_1 = list.collect{id, state -> state.fastqc_html_1}
        def fastqc_html_2 = list.collect{id, state -> state.fastqc_html_2}
        def fastqc_zip_1 = list.collect{id, state -> state.fastqc_zip_1}
        def fastqc_zip_2 = list.collect{id, state -> state.fastqc_zip_2}
        def trim_html_1 = list.collect{id, state -> state.trim_html_1}
        def trim_html_2 = list.collect{id, state -> state.trim_html_2}
        def trim_zip_1 = list.collect{id, state -> state.trim_zip_1}
        def trim_zip_2 = list.collect{id, state -> state.trim_zip_2}
        def trim_log_1 = list.collect{id, state -> state.trim_log_1}
        def trim_log_2 = list.collect{id, state -> state.trim_log_2}
        def fastq_1 = list.collect{id, state -> state.fastq_1}
        def fastq_2 = list.collect{id, state -> state.fastq_2}
        def sortmerna_log = list.collect{id, state -> state.sortmerna_log}
        def star_multiqc = list.collect{id, state -> state.star_multiqc}
        def salmon_quant_merged = list.collect{id, state -> state.salmon_quant_results}
        def genome_bam_sorted = list.collect{id, state -> state.genome_bam_sorted}
        def genome_bam_index = list.collect{id, state -> state.genome_bam_index}
        def genome_bam_stats = list.collect{id, state -> state.genome_bam_stats}
        def genome_bam_flagstat = list.collect{id, state -> state.genome_bam_flagstat}
        def genome_bam_idxstats = list.collect{id, state -> state.genome_bam_idxstats}
        def transcriptome_bam_sorted = list.collect{id, state -> state.transcriptome_bam_sorted}
        def transcriptome_bam_index = list.collect{id, state -> state.transcriptome_bam_index}
        def transcriptome_bam_stats = list.collect{id, state -> state.transcriptome_bam_stats}
        def transcriptome_bam_flagstat = list.collect{id, state -> state.transcriptome_bam_flagstat}
        def transcriptome_bam_idxstats = list.collect{id, state -> state.transcriptome_bam_idxstats}
        def markduplicates_multiqc = list.collect{id, state -> state.markduplicates_metrics}
        def featurecounts_multiqc = list.collect{id, state -> state.featurecounts_multiqc}
        def stringtie_transcript_gtf = list.collect{id, state -> state.stringtie_transcript_gtf}
        def stringtie_coverage_gtf = list.collect{id, state -> state.stringtie_coverage_gtf}
        def stringtie_abundance = list.collect{id, state -> state.stringtie_abundance}
        def stringtie_ballgown = list.collect{id, state -> state.stringtie_ballgown}
        def featurecounts = list.collect{id, state -> state.featurecounts}
        def featurecounts_summary = list.collect{id, state -> state.featurecounts_summary}
        def bedgraph_forward = list.collect{id, state -> state.bedgraph_forward}
        def bedgraph_reverse = list.collect{id, state -> state.bedgraph_reverse}
        def bigwig_forward = list.collect{id, state -> state.bigwig_forward}
        def bigwig_reverse = list.collect{id, state -> state.bigwig_reverse}
        def preseq_output = list.collect{id, state -> state.preseq_output}
        def qualimap_output_dir = list.collect{id, state -> state.qualimap_output_dir}
        def dupradar_output_dup_intercept_mqc = list.collect{id, state -> state.dupradar_output_dup_intercept_mqc}
        def bamstat_output = list.collect{id, state -> state.bamstat_output}
        // def inferexperiment_multiqc = list.collect{id, state -> state.inferexperiment_multiqc}
        def inner_dist_output_freq = list.collect{id, state -> state.inner_dist_output_freq}
        def junction_annotation_output_log = list.collect{id, state -> state.junction_annotation_output_log}
        def junction_saturation_output_plot_r = list.collect{id, state -> state.junction_saturation_output_plot_r}
        def read_distribution_output = list.collect{id, state -> state.read_distribution_output}
        def read_duplication_output_duplication_rate_mapping = list.collect{id, state -> state.read_duplication_output_duplication_rate_mapping}
        def tin_output_summary = list.collect{id, state -> state.tin_output_summary}
        ["merged", [
            fastqc_html_1: fastqc_html_1,
            fastqc_html_2: fastqc_html_2,
            fastqc_zip_1: fastqc_zip_1,
            fastqc_zip_2: fastqc_zip_2,
            fastqc_zip: fastqc_zip_1 + fastqc_zip_2,
            trim_html_1: trim_html_1,
            trim_html_2: trim_html_2,
            trim_zip_1: trim_zip_1, 
            trim_zip_2: trim_zip_2, 
            trim_zip: trim_zip_1 + trim_zip_2, 
            trim_log_1: trim_log_1, 
            trim_log_2: trim_log_2, 
            fastq_1: fastq_1,
            fastq_2: fastq_2,
            // sortmerna_multiqc: sortmerna_log,
            star_multiqc: star_multiqc, 
            salmon_multiqc: salmon_quant_merged,
            genome_bam_sorted: genome_bam_sorted,
            genome_bam_index: genome_bam_index, 
            samtools_stats: genome_bam_stats,
            samtools_flagstat: genome_bam_flagstat,
            samtools_dxstats: genome_bam_idxstats,
            transcriptome_bam_sorted: transcriptome_bam_sorted, 
            transcriptome_bam_index: transcriptome_bam_index, 
            transcriptome_bam_stats: transcriptome_bam_stats, 
            transcriptome_bam_flagstat: transcriptome_bam_flagstat, 
            transcriptome_bam_idxstats: transcriptome_bam_idxstats,
            markduplicates_multiqc: markduplicates_multiqc,
            featurecounts_multiqc: featurecounts_multiqc,
            stringtie_transcript_gtf: stringtie_transcript_gtf,
            stringtie_coverage_gtf: stringtie_coverage_gtf,
            stringtie_abundance: stringtie_abundance,
            stringtie_ballgown: stringtie_ballgown, 
            featurecounts: featurecounts,
            featurecounts_summary: featurecounts_summary, 
            bedgraph_forward: bedgraph_forward,
            bedgraph_reverse: bedgraph_reverse,
            bigwig_forward: bigwig_forward,
            bigwig_reverse: bigwig_reverse,
            preseq_output: preseq_output,
            qualimap_output_dir: qualimap_output_dir,
            dupradar_output_dup_intercept_mqc: dupradar_output_dup_intercept_mqc,
            bamstat_output: bamstat_output,
            inner_dist_output_freq: inner_dist_output_freq,
            junction_annotation_output_log: junction_annotation_output_log,
            junctionsaturation_multiqc: junction_saturation_output_plot_r,
            read_distribution_output: read_distribution_output,
            read_duplication_output_duplication_rate_mapping: read_duplication_output_duplication_rate_mapping,
            tin_output_summary: tin_output_summary, 
            salmon_quant_merged: salmon_quant_merged, 
            fasta: list[1][-1].fasta,
            gtf: list[1][-1].gtf,
            transcript_fasta: list[1][-1].transcript_fasta,
            gene_bed: list[1][-1].gene_bed,
            star_index: list[1][-1].star_index,
            salmon_index: list[1][-1].salmon_index,
            bbsplit_index: list[1][-1].bbxplit_index, 
            gtf_extra_attributes: list[1][-1].gtf_extra_attributes, 
            gtf_group_features: list[1][-1].gtf_group_features,
            pca_header_multiqc: list[1][-1].pca_header_multiqc, 
            clustering_header_multiqc: list[1][-1].clustering_header_multiqc,
            deseq2_vst: list[1][-1].deseq2_vst, 
            extra_deseq2_args: list[1][-1].extra_deseq2_args,
            extra_deseq2_args2: list[1][-1].extra_deseq2_args2,
            skip_deseq2_qc: list[1][-1].skip_deseq2_qc 
          ] 
        ]
    } 

    | salmon_quant_merge_counts.run (
        fromState: [ 
            "salmon_quant_results": "salmon_quant_merged", 
            "gtf": "gtf", 
            "gtf_extra_attributes": "gtf_extra_attributes", 
            "gtf_group_features": "gtf_group_features"
        ],
        toState: [
            "tpm_gene": "tpm_gene",
            "counts_gene": "counts_gene",
            "counts_gene_length_scaled": "counts_gene_length_scaled",
            "counts_gene_scaled": "counts_gene_scaled", 
            "tpm_transcript": "tpm_transcript", 
            "counts_transcript": "counts_transcript", 
            "salmon_merged_summarizedexperiment": "salmon_merged_summarizedexperiment"
        ]
    )

    | deseq2_qc.run (
        runIf: { id, state -> !state.skip_qc && !state.skip_deseq2_qc },
        fromState: [
            "counts": "counts_gene_length_scaled",
            "pca_header_multiqc": "pca_header_multiqc", 
            "clustering_header_multiqc": "clustering_header_multiqc",
            "deseq2_vst": "deseq2_vst", 
            "extra_deseq2_args": "extra_deseq2_args",
            "extra_deseq2_args2": "extra_deseq2_args2"
        ], 
        toState: [
            "deseq2_output": "deseq2_output", 
            "deseq2_pca_multiqc": "pca_multiqc", 
            "deseq2_dists_multiqc": "dists_multiqc"
        ]
    )

    | multiqc.run (
        fromState: [
          "multiqc_custom_config": "multiqc_custom_config", 
          "multiqc_title": "multiqc_title", 
          "multiqc_logo": "multiqc_logo",
          "multiqc_methods_description": "multiqc_methods_description",
        //   "workflow_summary": "workflow_summary", 
        //   "fail_trimming_multiqc": "fail_trimming_multiqc", 
        //   "fail_mapping_multiqc": "fail_mapping_multiqc", 
        //   "fail_strand_multiqc": "fail_strand_multiqc", 
          "fastqc_raw_multiqc": "fastqc_zip",
          "fastqc_trim_multiqc": "trim_zip",
          "trim_log_multiqc": "trim_log",
        //   "sortmerna_multiqc": "sortmerna_multiqc", 
          "star_multiqc": "star_multiqc", 
          "salmon_multiqc": "salmon_multiqc", 
          "samtools_stats": "genome_bam_stats", 
          "samtools_flagstat": "genome_bam_flagstat", 
          "samtools_idxstats": "genome_bam_idxstats", 
          "markduplicates_multiqc": "markduplicates_multiqc", 
          "featurecounts_multiqc": "featurecounts_multiqc", 
          "aligner_pca_multiqc": "deseq2_pca_multiqc", 
          "aligner_clustering_multiqc": "deseq2_dists_multiqc", 
          "preseq_multiqc": "preseq_output", 
          "qualimap_multiqc": "qualimap_output_dir", 
          "dupradar_multiqc": "dupradar_output_dup_intercept_mqc", 
          "bamstat_multiqc": "bamstat_output", 
        //   "inferexperiment_multiqc": "inferexperiment_multiqc", 
          "innerdistance_multiqc": "inner_dist_output_freq", 
          "junctionannotation_multiqc": "junction_annotation_output_log", 
          "junctionsaturation_multiqc": "junctionsaturation_multiqc",
          "readdistribution_multiqc": "read_distribution_output",
          "readduplication_multiqc": "read_duplication_output_duplication_rate_mapping", 
          "tin_multiqc": "tin_output_summary", 
          "deseq2_pca_multiqc": "deseq2_pca_multiqc", 
          "deseq2_dists_multiqc": "deseq2_pca_multiqc"
        ], 
        toState: [
          "multiqc_report": "report", 
          "multiqc_data": "data",
          "multiqc_plots": "plots", 
          "multiqc_versions": "versions"
        ]
    )
        
    // output_ch = qc_ch.concat(analysis_ch)
    // | map { id, state -> ["merged", state, join_id: ] }
    | map { id, state -> 
        [id, state + state.remove(skip_deseq2_qc) + state.remove(extra_deseq2_args2), state.remove(extra_deseq2_args) + state.remove(deseq2_vst) + state.remove(clustering_header_multiqc) + stat.remove(pca_header_multiqc) + state.remove(gtf_group_features) + state.remove(gtf_extra_attributes) + state.remove(markduplicates_multiqc) + state.remove(featurecounts_multiqc) + state.remove(salmon_multiqc) + state.remove(star_multiqc) + state.remove(samtools_stats) + state.remove(samtools_flagstat) + state.remove(samtools_idxstats) + state.remove(fastqc_zip) + state.remove(trim_zip)]}
    | setState (
      [
        "output_fasta": "fasta", 
        "output_gtf": "gtf", 
        "output_transcript_fasta": "transcript_fasta", 
        "output_gene_bed": "gene_bed", 
        "output_bbsplit_index": "bbsplit_index", 
        "output_star_index": "star_index", 
        "output_salmon_index": "salmon_index", 
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
        "transcriptome_bam_sorted": "transcriptome_bam_sorted", 
        "transcriptome_bam_index": "transcriptome_bam_index", 
        "transcriptome_bam_stats": "transcriptome_bam_stats", 
        "transcriptome_bam_flagstat": "transcriptome_bam_flagstat", 
        "transcriptome_bam_idxstats": "transcriptome_bam_idxstats",
        "salmon_quant_results": "salmon_quant_merged",
        "stringtie_transcript_gtf": "stringtie_transcript_gtf",
        "stringtie_coverage_gtf": "stringtie_coverage_gtf",
        "stringtie_abundance": "stringtie_abundance",
        "stringtie_ballgown": "stringtie_ballgown", 
        "featurecounts": "featurecounts",
        "featurecounts_summary": "featurecounts_summary", 
        "bedgraph_forward": "bedgraph_forward",
        "bedgraph_reverse": "bedgraph_reverse",
        "bigwig_forward": "bigwig_forward",
        "bigwig_reverse": "bigwig_reverse",
        "multiqc_report": "multiqc_report", 
        "multiqc_data": "multiqc_data"
      ]
    )

    | niceView()


  emit:
    output_ch
}

// ===============================
// === start of test workflows ===
// ===============================

// workflow test_wf {

//   // allow changing the resources_test dir
//   params.resources_test = params.rootDir + "/resources_test"

//   // or when running from s3: params.resources_test = "s3://openpipelines-data/"
//   testParams = [
//     param_list: [
//     ]
//   ]

//   output_ch =
//     channelFromParams(testParams, config)
//       | view { "Input: $it" }
//       | run_wf
//       | view { output ->
//         assert output.size() == 2 : "outputs should contain two elements; [id, file]"
//         // ...
//         "Output: $output"
//       }
//       | toSortedList()
//       | map { output_list ->
//         assert output_list.size() == 1 : "output channel should contain one event"
//         // ...
//       }
  
// }

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
        "  Provide the '--bam_csi_index' parameter to use a CSI instead of BAI index.\n\n" +
        "  Please see:\n" +
        "  https://github.com/nf-core/rnaseq/issues/744\n" +
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

//
// Create MultiQC tsv custom content from a list of values
//
def multiqcTsvFromList(tsv_data, header) {
  def tsv_string = ""
  if (tsv_data.size() > 0) {
    tsv_string += "${header.join('\t')}\n"
    tsv_string += tsv_data.join('\n')
  }
  return tsv_string
}