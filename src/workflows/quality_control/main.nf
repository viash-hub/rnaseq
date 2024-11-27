workflow run_wf {
    
    take: 
        input_ch

    main: 

        qc_ch = input_ch

            // temporary fix to force assignment when alignment in skipped
            | map {it} 

            // Feature biotype QC using featureCounts
            | map { id, state -> 
                def biotype_in_gtf = biotypeInGtf(state.gtf, state.biotype)
                def attribute_type = state.gencode ? "gene_type" : state.featurecounts_group_type
                def strand = (state.strandedness == "forward") ? 1 : ((state.strandedness == "reverse") ? 2 : 0)
                [ id, state + [biotype_in_gtf: biotype_in_gtf, attribute_type: attribute_type, strand: strand] ]
            }

            | featurecounts.run (
                runIf: { id, state -> !state.skip_qc && !state.skip_biotype_qc && state.biotype && state.biotype_in_gtf && !state.skip_align },
                fromState: [
                    "paired": "paired",
                    "strand": "strand", 
                    "annotation": "gtf", 
                    "input": "genome_bam", 
                    "attribute_type": "attribute_type",
                    "feature_type": "featurecounts_feature_type",
                    "count_read_pairs": "paired"
                ],
                toState: [
                    "featurecounts": "counts",
                    "featurecounts_summary": "summary"
                ],
                args: [
                    both_aligned: true,
                    same_strand: true
                ]
            )

        | multiqc_custom_biotype.run (
            runIf: { id, state -> !state.skip_qc && !state.skip_biotype_qc && state.biotype && state.featurecounts && !state.skip_align },
            fromState: [
                "id": "id",
                "biocounts": "featurecounts", 
                "biotypes_header": "biotypes_header"
            ],
            toState: [ 
                "featurecounts_multiqc": "featurecounts_multiqc", 
                "featurecounts_rrna_multiqc": "featurecounts_rrna_multiqc"
            ]
        )
              
        | preseq_lcextrap.run (
            runIf: { id, state -> !state.skip_qc && !state.skip_preseq && !state.skip_align },
            fromState: [
                "paired": "paired",
                "input": "genome_bam",
                "extra_preseq_args": "extra_preseq_args"
            ],
            toState: [ "preseq_output": "output" ]
        )
   
        | rseqc_bamstat.run (
            runIf: { id, state -> !state.skip_qc && !state.skip_rseqc && "bam_stat" in state.rseqc_modules && !state.skip_align },
            fromState: [
                "input_file": "genome_bam",
                "mapq": "map_qual"
            ],
            toState: [ "bamstat_output": "output" ]
        )
        | rseqc_inferexperiment.run(
            runIf: { id, state -> !state.skip_qc && !state.skip_rseqc && "infer_experiment" in state.rseqc_modules && !state.skip_align },
            fromState: [
                "input_file": "genome_bam",
                "refgene": "gene_bed",
                "sample_size": "sample_size",
                "mapq": "map_qual" 
            ],
            toState: [ "strandedness_output": "output" ]
        )
        // Get predicted strandedness from the RSeQC infer_experiment.py output
        | map { id, state -> 
            def inferred_strand = getInferexperimentStrandedness(state.strandedness_output, 30)
            def passed_strand_check = (state.strandedness != inferred_strand[0]) ? false : true
            [ id, state + [ inferred_strand: inferred_strand, passed_strand_check: passed_strand_check ] ]
        }
        | rseqc_inner_distance.run(
            runIf: { id, state -> !state.skip_qc && !state.skip_rseqc && state.paired && "inner_distance" in state.rseqc_modules && !state.skip_align },
            key: "inner_distance",
            fromState: [
                "input_file": "genome_bam",
                "refgene": "gene_bed",
                "sample_size": "sample_size",
                "mapq": "map_qual",
                "lower_bound": "lower_bound_size",
                "upper_bound": "upper_bound_size",
                "step": "step_size"
            ],
            toState: [ 
                "inner_dist_output_stats": "output_stats",
                "inner_dist_output_dist": "output_dist",
                "inner_dist_output_freq": "output_freq",
                "inner_dist_output_plot": "output_plot",
                "inner_dist_output_plot_r": "output_plot_r"
            ]
        )
        | rseqc_junctionannotation.run(
            runIf: { id, state -> !state.skip_qc && !state.skip_rseqc && "junction_annotation" in state.rseqc_modules && !state.skip_align },
            fromState: [
                "input": "genome_bam",
                "refgene": "gene_bed",
                "map_qual": "map_qual",
                "min_intron": "min_intron"
            ],
            toState: [
                "junction_annotation_output_log": "output_log",
                "junction_annotation_output_plot_r": "output_plot_r",
                "junction_annotation_output_junction_bed": "output_junction_bed",
                "junction_annotation_output_junction_interact": "output_junction_interact",
                "junction_annotation_output_junction_sheet": "output_junction_sheet",
                "junction_annotation_output_splice_events_plot": "output_splice_events_plot",
                "junction_annotation_output_splice_junctions_plot": "output_splice_junctions_plot" 
            ]
        )
        | rseqc_junctionsaturation.run(
            runIf: { id, state -> !state.skip_qc && !state.skip_rseqc && "junction_saturation" in state.rseqc_modules && !state.skip_align },
            fromState: [
                "input": "genome_bam",
                "refgene": "gene_bed",
                "sampling_percentile_lower_bound": "sampling_percentile_lower_bound",
                "sampling_percentile_upper_bound": "sampling_percentile_upper_bound",
                "sampling_percentile_step": "sampling_percentile_step",
                "min_intron": "min_intron",
                "min_splice_read": "min_splice_read",
                "map_qual": "map_qual"
            ],
            toState: [
                "junction_saturation_output_plot_r": "output_plot_r",
                "junction_saturation_output_plot": "output_plot"
            ]
        )
        | rseqc_readdistribution.run(
            runIf: { id, state -> !state.skip_qc && !state.skip_rseqc && "read_distribution" in state.rseqc_modules && !state.skip_align },
            fromState: [
                "input": "genome_bam",
                "refgene": "gene_bed", 
            ],
            toState: [ "read_distribution_output": "output" ]
        )                               
        | rseqc_readduplication.run(
            runIf: { id, state -> !state.skip_qc && !state.skip_rseqc && "read_duplication" in state.rseqc_modules && !state.skip_align },
            fromState: [
                "input": "genome_bam",
                "read_count_upper_limit": "read_count_upper_limit",
                "map_qual": "map_qual"
            ],
            toState: [
                "read_duplication_output_duplication_rate_plot_r": "output_duplication_rate_plot_r",
                "read_duplication_output_duplication_rate_plot": "output_duplication_rate_plot",
                "read_duplication_output_duplication_rate_mapping": "output_duplication_rate_mapping",
                "read_duplication_output_duplication_rate_sequence": "output_duplication_rate_sequence"
            ]
        )
        | rseqc_tin.run(
            runIf: { id, state -> !state.skip_qc && !state.skip_rseqc && "tin" in state.rseqc_modules && !state.skip_align },
            fromState: [
                "bam_input": "genome_bam",
                "bai_input": "genome_bam_index",
                "refgene": "gene_bed",
                "minimum_coverage": "minimum_coverage",
                "sample_size": "tin_sample_size",
                "subtract_background": "subtract_background"
            ],
            toState: [
                "tin_output_summary": "output_tin_summary",
                "tin_output_metrics": "output_tin"
            ]
        )

            | dupradar.run(
                runIf: { id, state -> !state.skip_qc && !state.skip_dupradar && !state.skip_align },
                fromState: [
                    "id": "id",
                    "input": "genome_bam",
                    "gtf_annotation": "gtf",
                    "paired": "paired",
                    "strandedness": "strandedness"
                ],
                toState: [ 
                    "dupradar_output_dupmatrix": "output_dupmatrix",
                    "dupradar_output_dup_intercept_mqc": "output_dup_intercept_mqc",
                    "dupradar_output_duprate_exp_boxplot": "output_duprate_exp_boxplot",
                    "dupradar_output_duprate_exp_densplot": "output_duprate_exp_densplot",
                    "dupradar_output_duprate_exp_denscurve_mqc": "output_duprate_exp_denscurve_mqc",
                    "dupradar_output_expression_histogram": "output_expression_histogram",
                    "dupradar_output_intercept_slope": "output_intercept_slope"
                ]
            )

        // TODO: Add outdir as an output argument to the qualimap module on biobox. 
        // Qualimap ouputs a few more raw data files to outdir but since the module is using a temporary directory as output dir these files are lost.
        | qualimap_rnaseq.run(
            fromState: [
                "bam": "genome_bam",
                "gtf": "gtf",
                "num_pr_bases": "pr_bases",
                "num_tr_bias": "tr_bias",
                "algorithm": "algorithm",
                "sequencing_protocol": "sequencing_protocol",
                "sorted": "sorted",
                "java_memory_size": "java_memory_size", 
            ],
            toState: [
                "qualimap_report": "report", 
                "qualimap_qc_report": "qc_report",
                "qualimap_counts": "counts"
            ]
        ) 

        merged_ch = qc_ch
            | toSortedList
            | map { list -> 
                def ids = list.collect { id, state -> id }
                def strandedness = list.collect { id, state -> state.strandedness }
                def num_trimmed_reads = list.collect { id, state -> state.num_trimmed_reads }
                def passed_trimmed_reads = list.collect { id, state -> state.passed_trimmed_reads }
                def passed_mapping = list.collect { id, state -> state.passed_mapping }
                def percent_mapped = list.collect { id, state -> state.percent_mapped }
                def inferred_strand = list.collect { id, state -> state.inferred_strand }
                def passed_strand_check = list.collect { id, state -> state.passed_strand_check }
                def gtf = list.collect { id, state -> state.gtf }.unique()[0]
                def gtf_extra_attributes = list.collect { id, state -> state.gtf_extra_attributes }.unique()[0]
                def gtf_group_features = list.collect { id, state -> state.gtf_group_features }.unique()[0]
                def pca_header_multiqc = list.collect { id, state -> state.pca_header_multiqc }.unique()[0]
                def clustering_header_multiqc = list.collect { id, state -> state.clustering_header_multiqc }.unique()[0]
                def aligner = list.collect { id, state -> state.aligner }.unique()[0]
                def pseudo_aligner = list.collect { id, state -> state.pseudo_aligner }.unique()[0]
                def deseq2_vst = list.collect { id, state -> state.deseq2_vst }.unique()[0]
                def extra_deseq2_args = list.collect { id, state -> state.extra_deseq2_args }.unique()[0]
                def extra_deseq2_args2 = list.collect { id, state -> state.extra_deseq2_args2 }.unique()[0]
                def skip_deseq2_qc = list.collect { id, state -> state.skip_deseq2_qc }.unique()[0] 
                def skip_qc = list.collect { id, state -> state.skip_qc }.unique()[0] 
                def skip_align = list.collect { id, state -> state.skip_align }.unique()[0] 
                def skip_pseudo_align = list.collect { id, state -> state.skip_pseudo_align }.unique()[0] 
                def quant_results = list.collect { id, state -> 
                    (state.quant_results_file instanceof java.nio.file.Path && state.quant_results_file.exists()) ? 
                        state.quant_results_file : 
                        null }
                def rsem_counts_gene = list.collect { id, state -> 
                    (state.rsem_counts_gene instanceof java.nio.file.Path && state.rsem_counts_gene.exists()) ? 
                        state.rsem_counts_gene : 
                        null }
                def rsem_counts_transcripts = list.collect { id, state -> 
                    (state.rsem_counts_transcripts instanceof java.nio.file.Path && state.rsem_counts_transcripts.exists()) ? 
                        state.rsem_counts_transcripts : 
                        null }
                def pseudo_quant_out_dir = list.collect { id, state -> 
                    (state.pseudo_quant_out_dir instanceof java.nio.file.Path && state.pseudo_quant_out_dir.exists()) ? 
                        state.pseudo_quant_out_dir : 
                        null }
                def pseudo_salmon_quant_results = list.collect { id, state -> 
                    (state.pseudo_salmon_quant_results_file instanceof java.nio.file.Path && state.pseudo_salmon_quant_results_file.exists()) ? 
                        state.pseudo_salmon_quant_results_file : 
                        null }
                def pseudo_kallisto_quant_results = list.collect { id, state -> 
                    (state.pseudo_kallisto_quant_results_file instanceof java.nio.file.Path && state.pseudo_kallisto_quant_results_file.exists()) ? 
                        state.pseudo_kallisto_quant_results_file : 
                        null }
                def fastqc_zip_1 = list.collect { id, state -> 
                    (state.fastqc_zip_1 instanceof java.nio.file.Path && state.fastqc_zip_1.exists()) ? 
                        state.fastqc_zip_1 : 
                        null }
                def fastqc_zip_2 = list.collect { id, state -> 
                    (state.fastqc_zip_2 instanceof java.nio.file.Path && state.fastqc_zip_2.exists()) ? 
                        state.fastqc_zip_2 : 
                        null }
                def trim_zip_1 = list.collect { id, state -> 
                    (state.trim_zip_1 instanceof java.nio.file.Path && state.trim_zip_1.exists()) ? 
                        state.trim_zip_1 : 
                        null }
                def trim_zip_2 = list.collect { id, state -> 
                    (state.trim_zip_2 instanceof java.nio.file.Path && state.trim_zip_2.exists()) ? 
                    state.trim_zip_2 : 
                    null }
                def trim_log_1 = list.collect { id, state -> 
                    (state.trim_log_1 instanceof java.nio.file.Path && state.trim_log_1.exists()) ? 
                    state.trim_log_1 : 
                    null }
                def trim_log_2 = list.collect { id, state -> 
                    (state.trim_log_2 instanceof java.nio.file.Path && state.trim_log_2.exists()) ? 
                        state.trim_log_2 : 
                        null }
                def sortmerna_multiqc = list.collect { id, state -> 
                    (state.sortmerna_multiqc instanceof java.nio.file.Path && state.sortmerna_multiqc.exists()) ? 
                        state.sortmerna_multiqc : 
                        null }
                def star_multiqc = list.collect { id, state -> 
                    (state.star_multiqc instanceof java.nio.file.Path && state.star_multiqc.exists()) ? 
                        state.star_multiqc : 
                        null }
                def genome_bam_stats = list.collect { id, state -> 
                    (state.genome_bam_stats instanceof java.nio.file.Path && state.genome_bam_stats.exists()) ? 
                        state.genome_bam_stats : 
                        null }
                def genome_bam_flagstat = list.collect { id, state -> 
                    (state.genome_bam_flagstat instanceof java.nio.file.Path && state.genome_bam_flagstat.exists()) ? 
                    state.genome_bam_flagstat : 
                    null }
                def genome_bam_idxstats = list.collect { id, state -> 
                    (state.genome_bam_idxstats instanceof java.nio.file.Path && state.genome_bam_idxstats.exists()) ? 
                        state.genome_bam_idxstats : 
                        null }
                def markduplicates_multiqc = list.collect { id, state -> 
                    (state.markduplicates_multiqc instanceof java.nio.file.Path && state.markduplicates_multiqc.exists()) ? 
                        state.markduplicates_multiqc : 
                        null }
                def salmon_multiqc = list.collect { id, state -> 
                    (state.salmon_multiqc instanceof java.nio.file.Path && state.salmon_multiqc.exists()) ? 
                        state.salmon_multiqc : 
                        null }
                def rsem_multiqc = list.collect { id, state -> 
                    (state.rsem_multiqc instanceof java.nio.file.Path && state.rsem_multiqc.exists()) ? 
                        state.rsem_multiqc : 
                        null }
                def pseudo_multiqc = list.collect { id, state -> 
                    (state.pseudo_multiqc instanceof java.nio.file.Path && state.pseudo_multiqc.exists()) ? 
                        state.pseudo_multiqc : 
                        null }
                def featurecounts_multiqc = list.collect { id, state -> 
                    (state.featurecounts_multiqc instanceof java.nio.file.Path && state.featurecounts_multiqc.exists()) ? 
                        state.featurecounts_multiqc : 
                        null }
                def featurecounts_rrna_multiqc = list.collect { id, state -> 
                    (state.featurecounts_rrna_multiqc instanceof java.nio.file.Path && state.featurecounts_rrna_multiqc.exists()) ? 
                        state.featurecounts_rrna_multiqc : 
                        null }
                def preseq_output = list.collect { id, state -> 
                    (state.preseq_output instanceof java.nio.file.Path && state.preseq_output.exists()) ? 
                        state.preseq_output : 
                        null }
                // def qualimap_output_dir = list.collect { id, state -> 
                //     (state.qualimap_output_dir instanceof java.nio.file.Path && state.qualimap_output_dir.exists()) ? 
                //         state.qualimap_output_dir : 
                //         null }
                def dupradar_output_dup_intercept_mqc = list.collect { id, state -> 
                    (state.dupradar_output_dup_intercept_mqc instanceof java.nio.file.Path && state.dupradar_output_dup_intercept_mqc.exists()) ? 
                        state.dupradar_output_dup_intercept_mqc : 
                            null }
                def dupradar_output_duprate_exp_denscurve_mqc = list.collect { id, state -> 
                    (state.dupradar_output_duprate_exp_denscurve_mqc instanceof java.nio.file.Path && state.dupradar_output_duprate_exp_denscurve_mqc.exists()) ? 
                        state.dupradar_output_duprate_exp_denscurve_mqc : 
                            null }
                def bamstat_output = list.collect { id, state -> 
                    (state.bamstat_output instanceof java.nio.file.Path && state.bamstat_output.exists()) ? 
                        state.bamstat_output : 
                        null }
                def inferexperiment_multiqc = list.collect { id, state -> 
                    (state.strandedness_output instanceof java.nio.file.Path && state.strandedness_output.exists()) ? 
                        state.strandedness_output : 
                        null }
                def inner_dist_output_freq = list.collect { id, state -> 
                    (state.inner_dist_output_freq instanceof java.nio.file.Path && state.inner_dist_output_freq.exists()) ? 
                        state.inner_dist_output_freq : 
                        null }
                def junction_annotation_output_log = list.collect { id, state -> 
                    (state.junction_annotation_output_log instanceof java.nio.file.Path && state.junction_annotation_output_log.exists()) ? 
                        state.junction_annotation_output_log : 
                        null }
                def junction_saturation_output_plot_r = list.collect { id, state -> 
                    (state.junction_saturation_output_plot_r instanceof java.nio.file.Path && state.junction_saturation_output_plot_r.exists()) ? 
                        state.junction_saturation_output_plot_r : 
                        null }
                def read_distribution_output = list.collect { id, state -> 
                    (state.read_distribution_output instanceof java.nio.file.Path && state.read_distribution_output.exists()) ? 
                        state.read_distribution_output : 
                        null }
                def read_duplication_output_duplication_rate_mapping = list.collect { id, state -> 
                    (state.read_duplication_output_duplication_rate_mapping instanceof java.nio.file.Path && state.read_duplication_output_duplication_rate_mapping.exists()) ? 
                        state.read_duplication_output_duplication_rate_mapping : 
                        null }
                def tin_output_summary = list.collect { id, state -> 
                    (state.tin_output_summary instanceof java.nio.file.Path && state.tin_output_summary.exists()) ? 
                        state.tin_output_summary : 
                        null }
                def multiqc_custom_config = list.collect { id, state -> state.multiqc_custom_config }.unique()[0] 
                ["merged", [
                    ids: ids, 
                    strandedness: strandedness, 
                    num_trimmed_reads: num_trimmed_reads,
                    passed_trimmed_reads: passed_trimmed_reads,
                    passed_mapping: passed_mapping,
                    percent_mapped: percent_mapped,
                    inferred_strand: inferred_strand, 
                    passed_strand_check: passed_strand_check, 
                    skip_align: skip_align,
                    skip_pseudo_align: skip_pseudo_align,
                    quant_results: quant_results.findAll { it != null }, 
                    rsem_counts_gene: rsem_counts_gene.findAll { it != null },
                    rsem_counts_transcripts: rsem_counts_transcripts.findAll { it != null },
                    pseudo_quant_out_dir: pseudo_quant_out_dir.findAll { it != null },
                    pseudo_salmon_quant_results: pseudo_salmon_quant_results.findAll { it != null },
                    pseudo_kallisto_quant_results: pseudo_kallisto_quant_results.findAll { it != null },
                    gtf: gtf, 
                    gtf_extra_attributes: gtf_extra_attributes, 
                    gtf_group_features: gtf_group_features,
                    pca_header_multiqc: pca_header_multiqc, 
                    clustering_header_multiqc: clustering_header_multiqc,
                    aligner: aligner,
                    pseudo_aligner: pseudo_aligner,
                    deseq2_vst: deseq2_vst, 
                    extra_deseq2_args: extra_deseq2_args,
                    extra_deseq2_args2: extra_deseq2_args2,
                    skip_deseq2_qc: skip_deseq2_qc,
                    fastqc_zip: fastqc_zip_1 + fastqc_zip_2,
                    trim_zip: trim_zip_1 + trim_zip_2, 
                    trim_log: trim_log_1 + trim_log_2, 
                    sortmerna_multiqc: sortmerna_multiqc,
                    star_multiqc: star_multiqc, 
                    genome_bam_stats: genome_bam_stats,
                    genome_bam_flagstat: genome_bam_flagstat,
                    genome_bam_idxstats: genome_bam_idxstats,
                    markduplicates_multiqc: markduplicates_multiqc,
                    salmon_multiqc: salmon_multiqc,
                    rsem_multiqc: rsem_multiqc,
                    pseudo_multiqc: pseudo_multiqc,
                    featurecounts_multiqc: featurecounts_multiqc,
                    featurecounts_rrna_multiqc: featurecounts_rrna_multiqc,
                    preseq_output: preseq_output,
                    // qualimap_output_dir: qualimap_output_dir,
                    dupradar_output_dup_intercept_mqc: dupradar_output_dup_intercept_mqc,
                    dupradar_output_duprate_exp_denscurve_mqc: dupradar_output_duprate_exp_denscurve_mqc,
                    bamstat_output: bamstat_output,
                    inner_dist_output_freq: inner_dist_output_freq,
                    inferexperiment_multiqc: inferexperiment_multiqc,
                    junction_annotation_output_log: junction_annotation_output_log,
                    junction_saturation_output_plot_r: junction_saturation_output_plot_r,
                    read_distribution_output: read_distribution_output,
                    read_duplication_output_duplication_rate_mapping: read_duplication_output_duplication_rate_mapping,
                    tin_output_summary: tin_output_summary, 
                    multiqc_custom_config: multiqc_custom_config
                ] ]
            } 

            // Merge quantification results of alignment
            | merge_quant_results.run (
                runIf: { id, state -> !state.skip_align && state.aligner == 'star_salmon' },
                fromState: [ 
                    "salmon_quant_results": "quant_results", 
                    "gtf": "gtf", 
                    "gtf_extra_attributes": "gtf_extra_attributes", 
                    "gtf_group_features": "gtf_group_features"
                ],
                args: [ quant_type: "salmon"],
                toState: [
                    "tpm_gene": "tpm_gene",
                    "counts_gene": "counts_gene",
                    "counts_gene_length_scaled": "counts_gene_length_scaled",
                    "counts_gene_scaled": "counts_gene_scaled", 
                    "tpm_transcript": "tpm_transcript", 
                    "counts_transcript": "counts_transcript", 
                    "lengths_gene": "lengths_gene",
                    "lengths_transcript": "lengths_transcript",
                    "quant_merged_summarizedexperiment": "quant_merged_summarizedexperiment"
                ], 
                key: "merge_quant_results"
            )

            | rsem_merge_counts.run (
                runIf: { id, state -> state.aligner == 'star_rsem' }, 
                fromState: [
                    "counts_gene": "rsem_counts_gene",
                    "counts_transcripts": "rsem_counts_transcripts"
                ],
                toState: [
                    "tpm_gene": "merged_gene_tpm",
                    "counts_gene": "merged_gene_counts",
                    "tpm_transcript": "merged_transcript_tpm", 
                    "counts_transcript": "merged_transcript_counts"
                ]
            )

            | deseq2_qc.run (
                runIf: { id, state -> !state.skip_qc && !state.skip_deseq2_qc && !state.skip_align },
                fromState: { id, state ->
                    def counts = (state.aligner == "star_rsem") ? state.counts_gene : state.counts_gene_length_scaled
                    [
                        counts: counts,
                        vst: state.deseq2_vst, 
                        label: state.aligner 
                    ]
                },
                args: [count_col: 3, id_col: 1, outprefix: "deseq2"], 
                toState: [
                    "deseq2_output": "outdir", 
                    "deseq2_pca_multiqc": "pca_multiqc", 
                    "deseq2_dists_multiqc": "dists_multiqc" 
                ], 
                key: "deseq2_qc_align_quant"
            )

            // Merge quantification results of pseudo alignment
            | merge_quant_results.run (
                runIf: { id, state -> !state.skip_pseudo_align },
                fromState: [ 
                    "salmon_quant_results": "pseudo_salmon_quant_results",
                    "kallisto_quant_results": "pseudo_kallisto_quant_results",
                    "gtf": "gtf", 
                    "gtf_extra_attributes": "gtf_extra_attributes", 
                    "gtf_group_features": "gtf_group_features",
                    "quant_type": "pseudo_aligner"
                ],
                toState: [
                    "pseudo_tpm_gene": "tpm_gene",
                    "pseudo_counts_gene": "counts_gene",
                    "pseudo_counts_gene_length_scaled": "counts_gene_length_scaled",
                    "pseudo_counts_gene_scaled": "counts_gene_scaled", 
                    "pseudo_tpm_transcript": "tpm_transcript", 
                    "pseudo_counts_transcript": "counts_transcript", 
                    "pseudo_lengths_gene": "lengths_gene",
                    "pseudo_lengths_transcript": "lengths_transcript",
                    "pseudo_quant_merged_summarizedexperiment": "quant_merged_summarizedexperiment"
                ], 
                key: "merge_pseudo_quant_results"
            )

            | deseq2_qc.run (
                runIf: { id, state -> !state.skip_qc && !state.skip_deseq2_qc && !state.skip_pseudo_align },
                fromState: [
                    "counts": "pseudo_counts_gene_length_scaled",
                    "vst": "deseq2_vst", 
                    "label": "pseudo_aligner" 
                ],
                args: [count_col: 3, id_col: 1, outprefix: "deseq2"], 
                toState: [
                    "deseq2_output_pseudo": "outdir", 
                    "deseq2_pca_multiqc_pseudo": "pca_multiqc", 
                    "deseq2_dists_multiqc_pseudo": "dists_multiqc" 
                ],
                key: "deseq2_qc_pseuso_align_quant"
            )
            
            // Get list of samples that failed trimming, mapping, and strand check for MultiQC report
            | map { id, state -> 
                def fail_trimming_header = ["Sample", "Reads after trimming"]
                def fail_trimming_multiqc = ""
                def star_mapping_header = ["Sample", "STAR uniquely mapped reads (%)"] 
                def fail_mapping_multiqc = ""
                def strand_check_header = ["Sample", "Provided strandedness", "Inferred strandedness", "Sense (%)", "Antisense (%)", "Undetermined (%)"]
                def fail_strand_multiqc = ""
                if (state.ids.size() > 0) {
                    fail_trimming_multiqc += "${fail_trimming_header.join('\t')}\n"
                    fail_mapping_multiqc += "${star_mapping_header.join('\t')}\n"
                    fail_strand_multiqc += "${strand_check_header.join('\t')}\n"
                    for (i=0; i<state.ids.size(); i++) {
                        if (!state.passed_trimmed_reads[i]) {
                            tsv_data = [state.ids[i], state.num_trimmed_reads[i]].join('\t')
                            fail_trimming_multiqc += tsv_data.join('\n')
                        }
                        if (!state.passed_mapping[i]) {
                            tsv_data = [state.ids[i], state.percent_mapped[i]].join('\t')
                            fail_mapping_multiqc += tsv_data.join('\n')
                        }
                        if (!state.passed_strand_check[i]) {
                            tsv_data = ([state.ids[i], state.strandedness[i]] + state.inferred_strand[i]).join('\t')
                            fail_strand_multiqc += tsv_data.join('\n')
                        }
                    }
                }
                
                [ id, state + [fail_trimming_multiqc: fail_trimming_multiqc, fail_mapping_multiqc: fail_mapping_multiqc, fail_strand_multiqc: fail_strand_multiqc] ]
            }

            | map { id, state -> 
                state.each { key, value ->
                    if (value instanceof ArrayList) {
                        value.removeAll { it == null }
                    }
                }
                mod_state = state.findAll { key, value -> value != null }
                [ id, mod_state ]
            }

            | prepare_multiqc_input.run(
                runIf: { id, state -> !state.skip_qc && !state.skip_multiqc },
                fromState: [
                    "fail_trimming_multiqc": "fail_trimming_multiqc", 
                    "fail_mapping_multiqc": "fail_mapping_multiqc", 
                    "fail_strand_multiqc": "fail_strand_multiqc", 
                    "fastqc_raw_multiqc": "fastqc_zip",
                    "fastqc_trim_multiqc": "trim_zip",
                    "trim_log_multiqc": "trim_log",
                    "sortmerna_multiqc": "sortmerna_multiqc", 
                    "star_multiqc": "star_multiqc", 
                    "salmon_multiqc": "salmon_multiqc", 
                    "rsem_multiqc": "rsem_multiqc", 
                    "pseudo_multiqc": "pseudo_multiqc",
                    "samtools_stats": "genome_bam_stats", 
                    "samtools_flagstat": "genome_bam_flagstat", 
                    "samtools_idxstats": "genome_bam_idxstats", 
                    "markduplicates_multiqc": "markduplicates_multiqc",
                    "featurecounts_multiqc": "featurecounts_multiqc",
                    "featurecounts_rrna_multiqc": "featurecounts_rrna_multiqc", 
                    "aligner_pca_multiqc": "deseq2_pca_multiqc", 
                    "aligner_clustering_multiqc": "deseq2_dists_multiqc", 
                    "pseudo_aligner_pca_multiqc": "deseq2_pca_multiqc_pseudo", 
                    "pseudo_aligner_clustering_multiqc": "deseq2_dists_multiqc_pseudo", 
                    "preseq_multiqc": "preseq_output", 
                    // "qualimap_multiqc": "qualimap_output_dir", 
                    "dupradar_output_dup_intercept_mqc": "dupradar_output_dup_intercept_mqc", 
                    "dupradar_output_duprate_exp_denscurve_mqc": "dupradar_output_duprate_exp_denscurve_mqc",
                    "bamstat_multiqc": "bamstat_output", 
                    "inferexperiment_multiqc": "inferexperiment_multiqc", 
                    "innerdistance_multiqc": "inner_dist_output_freq", 
                    "junctionannotation_multiqc": "junction_annotation_output_log", 
                    "junctionsaturation_multiqc": "junction_saturation_output_plot_r",
                    "readdistribution_multiqc": "read_distribution_output",
                    "readduplication_multiqc": "read_duplication_output_duplication_rate_mapping", 
                    "tin_multiqc": "tin_output_summary",
                    "multiqc_config": "multiqc_custom_config"
                ], 
                toState: [ "multiqc_input": "output" ]
            )

            | multiqc.run (
                runIf: { id, state -> !state.skip_qc && !state.skip_multiqc },
                fromState: [
                    "title": "multiqc_title", 
                    "input": "multiqc_input", 
                ], 
                args: [exclude_modules: "general_stats"],
                toState: [
                    "multiqc_report": "output_report", 
                    "multiqc_data": "output_data",
                    "multiqc_plots": "output_plots"
                ]
            )

            | map { id, state -> 
                [ id, [ 
                    tpm_gene: state.tpm_gene,
                    counts_gene: state.counts_gene,
                    counts_gene_length_scaled: state.counts_gene_length_scaled,
                    counts_gene_scaled: state.counts_gene_scaled, 
                    tpm_transcript: state.tpm_transcript, 
                    counts_transcript: state.counts_transcript, 
                    quant_merged_summarizedexperiment: state.quant_merged_summarizedexperiment,
                    deseq2_output: state.deseq2_output, 
                    pseudo_tpm_gene: state.pseudo_tpm_gene,
                    pseudo_counts_gene: state.pseudo_counts_gene,
                    pseudo_counts_gene_length_scaled: state.pseudo_counts_gene_length_scaled,
                    pseudo_counts_gene_scaled: state.pseudo_counts_gene_scaled, 
                    pseudo_tpm_transcript: state.pseudo_tpm_transcript, 
                    pseudo_counts_transcript: state.pseudo_counts_transcript, 
                    pseudo_quant_merged_summarizedexperiment: state.pseudo_quant_merged_summarizedexperiment,
                    deseq2_output_pseudo: state.deseq2_output_pseudo,
                    multiqc_report: state.multiqc_report, 
                    multiqc_data: state.multiqc_data, 
                    multiqc_plots: state.multiqc_plots
                ] ] 
            }

            | map { list -> list[1]}

        output_ch = qc_ch
        
            | combine(merged_ch)

            | map { list -> [list[0], list[1] + list[2]] }

            | map { id, state -> 
                def mod_state = state.findAll { key, value -> value instanceof java.nio.file.Path && value.exists() }
                [ id, mod_state ]
            }

            | setState (
                [
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
                    "qualimap_report": "qualimap_report", 
                    "qualimap_qc_report": "qualimap_qc_report",
                    "qualimap_counts": "qualimap_counts",
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
                    "lengths_gene": "lengths_gene",
                    "lengths_transcript": "lengths_transcript",
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

    emit:
        output_ch
}

//
// Function to check whether biotype field exists in GTF file
//
def biotypeInGtf(gtf_file, biotype) {
    def hits = 0
    gtf_file.eachLine { line ->
        def attributes = line.split('\t')[-1].split()
        if (attributes.contains(biotype)) {
            hits += 1
        }
    }
    if (hits) {
        return true
    } else {
        log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Biotype attribute '${biotype}' not found in the last column of the GTF file!\n\n" +
            "  Biotype QC will be skipped to circumvent the issue below:\n" +
            "  https://github.com/nf-core/rnaseq/issues/460\n\n" +
            "  Amend '--featurecounts_group_type' to change this behaviour.\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        return false
    }
}

//
// Function that parses and returns the predicted strandedness from the RSeQC infer_experiment.py output
//
def getInferexperimentStrandedness(inferexperiment_file, cutoff=30) {
    def sense        = 0
    def antisense    = 0
    def undetermined = 0
    inferexperiment_file.eachLine { line ->
        def undetermined_matcher = line =~ /Fraction of reads failed to determine:\s([\d\.]+)/
        def se_sense_matcher     = line =~ /Fraction of reads explained by "\++,--":\s([\d\.]+)/
        def se_antisense_matcher = line =~ /Fraction of reads explained by "\+-,-\+":\s([\d\.]+)/
        def pe_sense_matcher     = line =~ /Fraction of reads explained by "1\++,1--,2\+-,2-\+":\s([\d\.]+)/
        def pe_antisense_matcher = line =~ /Fraction of reads explained by "1\+-,1-\+,2\+\+,2--":\s([\d\.]+)/
        if (undetermined_matcher) undetermined = undetermined_matcher[0][1].toFloat() * 100
        if (se_sense_matcher)     sense        = se_sense_matcher[0][1].toFloat() * 100
        if (se_antisense_matcher) antisense    = se_antisense_matcher[0][1].toFloat() * 100
        if (pe_sense_matcher)     sense        = pe_sense_matcher[0][1].toFloat() * 100
        if (pe_antisense_matcher) antisense    = pe_antisense_matcher[0][1].toFloat() * 100
    }
    def strandedness = 'unstranded'
    if (sense >= 100-cutoff) {
        strandedness = 'forward'
    } else if (antisense >= 100-cutoff) {
        strandedness = 'reverse'
    }
    return [ strandedness, sense, antisense, undetermined ]
}



