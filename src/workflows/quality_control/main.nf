workflow run_wf {
    
    take: 
        input_ch

    main: 

        qc_ch = input_ch

        // Feature biotype QC using featureCounts
        | map { id, state -> 
            def biotype_in_gtf = biotypeInGtf(state.gtf, state.biotype)
            [ id, state + [biotype_in_gtf: biotype_in_gtf] ]
        }

        | subread_featurecounts.run (
            runIf: { id, state -> !state.skip_qc && !state.skip_biotype_qc && state.biotype && state.biotype_in_gtf },
            fromState: [
                "paired": "paired", 
                "strandedness": "strandedness", 
                "gtf": "gtf", 
                "bam": "genome_bam", 
                "gencode": "gencode",
                "extra_featurecounts_args": "extra_featurecounts_args",
                "featurecounts_group_type": "featurecounts_group_type",
                "featurecounts_feature_type": "featurecounts_feature_type",
            ],
            toState: [
                "featurecounts": "counts",
                "featurecounts_summary": "summary"
            ]
        )

        | multiqc_custom_biotype.run (
            runIf: { id, state -> !state.skip_qc && !state.skip_biotype_qc && state.biotype && state.featurecounts },
            fromState: [
                "id": "id",
                "biocounts": "featurecounts", 
                "biotypes_header": "biotypes_header"
            ],
            toState: [ "featurecounts_multiqc": "featurecounts_multiqc" ]
        )
              
        | preseq_lcextrap.run (
            runIf: { id, state -> !state.skip_qc && !state.skip_preseq },
            fromState: [
                "paired": "paired",
                "bam": "genome_bam",
                "extra_preseq_args": "extra_preseq_args"
            ],
            toState: [ "preseq_output": "output" ],
        )
   
        | rseqc_bamstat.run (
            runIf: { id, state -> "bam_stat" in state.rseqc_modules },
            fromState: [
                "input": "genome_bam",
                "map_qual": "map_qual"
            ],
            toState: [ "bamstat_output": "output" ]
        )
        | rseqc_inferexperiment.run(
            runIf: { id, state -> "infer_experiment" in state.rseqc_modules },
            fromState: [
                "input": "genome_bam",
                "refgene": "gene_bed",
                "sample_size": "sample_size",
                "map_qual": "map_qual"
            ],
            toState: [ "strandedness_output": "output" ]
        )
        // TODO: Get predicted strandedness from the RSeQC infer_experiment.py output
        | rseqc_innerdistance.run(
            runIf: { id, state -> state.paired && "inner_distance" in state.rseqc_modules },
            key: "inner_distance",
            fromState: [
                "input": "genome_bam",
                "refgene": "gene_bed",
                "paired": "paired",
                "sample_size": "sample_size",
                "map_qual": "map_qual",
                "lower_bound_size": "lower_bound_size",
                "upper_bound_size": "upper_bound_size",
                "step_size": "step_size"
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
            runIf: { id, state -> "junction_annotation" in state.rseqc_modules },
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
            runIf: { id, state -> "junction_saturation" in state.rseqc_modules },
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
            runIf: { id, state -> "read_distribution" in state.rseqc_modules },
            fromState: [
                "input": "genome_bam",
                "refgene": "gene_bed"
            ],
            toState: [ "read_distribution_output": "output" ]
        )                               
        | rseqc_readduplication.run(
            runIf: { id, state -> "read_duplication" in state.rseqc_modules },
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
            runIf: { id, state -> "tin" in state.rseqc_modules },
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

        | qualimap.run(
            fromState: [
                "input": "genome_bam",
                "gtf": "gtf",
                "pr_bases": "pr_bases",
                "tr_bias": "tr_bias",
                "algorithm": "algorithm",
                "sequencing_protocol": "sequencing_protocol",
                "sorted": "sorted",
                "java_memory_size": "java_memory_size"
            ],
            toState: [
                "qualimap_output_pdf": "output_pdf",
                "qualimap_output_dir": "output_dir"
            ]
        ) 

        merged_ch = qc_ch
        
        // TODO:
        // def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)
        // workflow_summary = WorkflowRnaseq.paramsSummaryMultiqc(workflow, summary_params)
        // methods_description = WorkflowRnaseq.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
        | toSortedList
        | map { list -> 
            def fastqc_zip_1 = list.collect{id, state -> state.fastqc_zip_1}
            def fastqc_zip_2 = list.collect{id, state -> state.fastqc_zip_2}
            def trim_zip_1 = list.collect{id, state -> state.trim_zip_1}
            def trim_zip_2 = list.collect{id, state -> state.trim_zip_2}
            def trim_log_1 = list.collect{id, state -> state.trim_log_1}
            def trim_log_2 = list.collect{id, state -> state.trim_log_2}
            def sortmerna_log = list.collect{id, state -> state.sortmerna_log}
            def star_multiqc = list.collect{id, state -> state.star_multiqc}
            def salmon_quant_results = list.collect{id, state -> state.salmon_quant_results}
            def genome_bam_stats = list.collect{id, state -> state.genome_bam_stats}
            def genome_bam_flagstat = list.collect{id, state -> state.genome_bam_flagstat}
            def genome_bam_idxstats = list.collect{id, state -> state.genome_bam_idxstats}
            def markduplicates_multiqc = list.collect{id, state -> state.markduplicates_multiqc}
            def featurecounts_multiqc = list.collect{id, state -> state.featurecounts_multiqc}
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
            def gtf = list.collect{id, state -> state.gtf}.unique()[0]
            def gtf_extra_attributes = list.collect{id, state -> state.gtf_extra_attributes}.unique()[0]
            def gtf_group_features = list.collect{id, state -> state.gtf_group_features}.unique()[0]
            def pca_header_multiqc = list.collect{id, state -> state.pca_header_multiqc}.unique()[0]
            def clustering_header_multiqc = list.collect{id, state -> state.clustering_header_multiqc}.unique()[0]
            def deseq2_vst = list.collect{id, state -> state}.deseq2_vst.unique()[0]
            def extra_deseq2_args = list.collect{id, state -> state.extra_deseq2_args}.unique()[0]
            def extra_deseq2_args2 = list.collect{id, state -> state.extra_deseq2_args2}.unique()[0]
            def skip_deseq2_qc = list.collect{id, state -> state.skip_deseq2_qc}.unique()[0] 
            ["merged", [
                fastqc_zip: fastqc_zip_1 + fastqc_zip_2,
                trim_zip: trim_zip_1 + trim_zip_2, 
                trim_log: trim_log_1 + trim_log_2, 
                sortmerna_log: sortmerna_log,
                star_multiqc: star_multiqc, 
                salmon_multiqc: salmon_quant_results,
                genome_bam_stats: genome_bam_stats,
                genome_bam_flagstat: genome_bam_flagstat,
                genome_bam_idxstats: genome_bam_idxstats,
                markduplicates_multiqc: markduplicates_multiqc,
                featurecounts_multiqc: featurecounts_multiqc,
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
                salmon_quant_results: salmon_quant_results, 
                gtf: gtf, 
                gtf_extra_attributes: gtf_extra_attributes, 
                gtf_group_features: gtf_group_features,
                pca_header_multiqc: pca_header_multiqc, 
                clustering_header_multiqc: clustering_header_multiqc,
                deseq2_vst: deseq2_vst, 
                extra_deseq2_args: extra_deseq2_args,
                extra_deseq2_args2: extra_deseq2_args2,
                skip_deseq2_qc: skip_deseq2_qc 
            ] ]
        } 

        | salmon_quant_merge_counts.run (
            fromState: [ 
                "salmon_quant_results": "salmon_quant_results", 
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
            ], 
            auto: [publish: true]
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
            ], 
            auto: [publish: true]
        )

        | prepare_multiqc_input.run(
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
            //   "sortmerna_multiqc": "sortmerna_log", 
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
              "multiqc_input": "output"
            ]
        )
        
        | multiqc.run (
            fromState: [
              "multiqc_custom_config": "multiqc_custom_config", 
              "multiqc_title": "multiqc_title", 
              "multiqc_logo": "multiqc_logo",
              "multiqc_methods_description": "multiqc_methods_description",
            //   "workflow_summary": "workflow_summary", 
              "input": "multiqc_input"
            ], 
            toState: [
              "multiqc_report": "report", 
              "multiqc_data": "data",
              "multiqc_plots": "plots", 
              "multiqc_versions": "versions"
            ], 
            auto: [publish: true]
        )
        
        output_ch = qc_ch
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
                "qualimap_output_dir": "qualimap_output_dir",
                "qualimap_output_pdf": "qualimap_output_pdf", 
                "featurecounts": "featurecounts",
                "featurecounts_summary": "featurecounts_summary" 
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