workflow run_wf {
    
    take: 
        input_ch

    main: 
        output_ch = input_ch 
  
        | 
            (        
                rseqc_bamstat.run (
                    fromState: [
                            "input": "bam_input",
                            "map_qual": "map_qual",
                            "output": "bamstat_output"
                            ],
                    toState: { id, output, state ->
                        [ "bamstat_output": output.output ]
                        }

                )

                & 
                
                rseqc_inferexperiment.run(
                    fromState: [
                            "input": "bam_input",
                            "refgene": "refgene",
                            "sample_size": "sample_size",
                            "map_qual": "map_qual",
                            "output": "strandedness_output"
                        ],
                    toState: { id, output, state ->
                        [ "strandedness_output": output.output ]
                        }
                )

                & 
                
                rseqc_innerdistance.run(
                    runIf: {id, state -> state.paired},
                    fromState: [
                            "input": "bam_input",
                            "refgene": "refgene",
                            "paired": "paired",
                            "sample_size": "sample_size",
                            "map_qual": "map_qual",
                            "lower_bound_size": "lower_bound_size",
                            "upper_bound_size": "upper_bound_size",
                            "step_size": "step_size",
                            "output_stats": "inner_dist_output_stats",
                            "output_dist": "inner_dist_output_dist",
                            "output_freq": "inner_dist_output_freq",
                            "output_plot": "inner_dist_output_plot",
                            "output_plot_r": "inner_dist_output_plot_r"
                        ],
                    toState: { id, output, state -> 
                        [
                            "inner_dist_output_stats": output.output_stats,
                            "inner_dist_output_dist": output.output_dist,
                            "inner_dist_output_freq": output.output_freq,
                            "inner_dist_output_plot": output.output_plot,
                            "inner_dist_output_plot_r": output.output_plot_r
                        ]}
                )

                & 
                
                rseqc_junctionannotation.run(
                    fromState: [
                            "input": "bam_input",
                            "refgene": "refgene",
                            "map_qual": "map_qual",
                            "min_intron": "min_intron",
                            "output_log": "junction_annotation_output_log",
                            "output_plot_r": "junction_annotation_output_plot_r",
                            "output_junction_bed": "junction_annotation_output_junction_bed",
                            "output_junction_interact": "junction_annotation_output_junction_interact",
                            "output_junction_sheet": "junction_annotation_output_junction_sheet",
                            "output_splice_events_plot": "junction_annotation_output_splice_events_plot",
                            "output_splice_junctions_plot": "junction_annotation_output_splice_junctions_plot"
                        ],
                    toState: { id, output, state ->
                        [
                            "junction_annotation_output_log": output.output_log,
                            "junction_annotation_output_plot_r": output.output_plot_r,
                            "junction_annotation_output_junction_bed": output.output_junction_bed,
                            "junction_annotation_output_junction_interact": output.output_junction_interact,
                            "junction_annotation_output_junction_sheet": output.output_junction_sheet,
                            "junction_annotation_output_splice_events_plot": output.output_splice_events_plot,
                            "junction_annotation_output_splice_junctions_plot": output.output_splice_junctions_plot
                        ]

                    }
                )

                & 

                rseqc_junctionsaturation.run(
                    fromState: [
                        "input": "bam_input",
                        "refgene": "refgene",
                        "sampling_percentile_lower_bound": "sampling_percentile_lower_bound",
                        "sampling_percentile_upper_bound": "sampling_percentile_upper_bound",
                        "sampling_percentile_step": "sampling_percentile_step",
                        "min_intron": "min_intron",
                        "min_splice_read": "min_splice_read",
                        "map_qual": "map_qual",
                        "output_plot_r": "junction_saturation_output_plot_r",
                        "output_plot": "junction_saturation_output_plot"

                    ],
                    toState: { id, output, state ->
                        [
                            "junction_saturation_output_plot_r": output.output_plot_r,
                            "junction_saturation_output_plot": output.output_plot
                        ]
                    }
                )

                &

                rseqc_readdistribution.run(
                    fromState: [
                            "input": "bam_input",
                            "refgene": "refgene",
                            "output": "read_distribution_output"
                        ],
                    toState: { id, output, state ->
                        [ "read_distribution_output": output.output ]
                        }
                )
                                
                &

                rseqc_readduplication.run(
                    fromState: [
                            "input": "bam_input",
                            "read_count_upper_limit": "read_count_upper_limit",
                            "map_qual": "map_qual",
                            "output_duplication_rate_plot_r": "read_duplication_output_duplication_rate_plot_r",
                            "output_duplication_rate_plot": "read_duplication_output_duplication_rate_plot",
                            "output_duplication_rate_mapping": "read_duplication_output_duplication_rate_mapping",
                            "output_duplication_rate_sequence": "read_duplication_output_duplication_rate_sequence"
                        ],
                    toState: { id, output, state ->
                        [
                             "read_duplication_output_duplication_rate_plot_r": output.output_duplication_rate_plot_r,
                             "output_duplication_rate_plot": output.read_duplication_output_duplication_rate_plot,
                             "output_duplication_rate_mapping": output.read_duplication_output_duplication_rate_mapping,
                             "output_duplication_rate_sequence": output.read_duplication_output_duplication_rate_sequence     
                        ]
                    }
                )

                &

                rseqc_tin.run(
                    fromState: [
                            "bam_input": "bam_input",
                            "bai_input": "bai_input",
                            "refgene": "refgene",
                            "minimum_coverage": "minimum_coverage",
                            "sample_size": "tin_sample_size",
                            "subtract_background": "subtract_background",
                            "output_tin_summary": "tin_output_summary",
                            "output_tin": "tin_output_metrics"
                        ],
                    toState: { id, output, state ->
                        [
                            "tin_output_summary": output.output_tin_summary,
                            "tin_output_metrics": output.output_tin,
                        ]
                    }
                )

                &

                dupradar.run(
                    fromState: [
                            "input": "bam_input",
                            "gtf_annotation": "gtf_input",
                            "paired": "paired",
                            "strandedness": "strandedness",
                            "output_dupmatrix": "dupradar_output_dupmatrix",
                            "output_dup_intercept_mqc": "dupradar_output_dup_intercept_mqc",
                            "output_duprate_exp_boxplot": "dupradar_output_duprate_exp_boxplot",
                            "output_duprate_exp_densplot": "dupradar_output_duprate_exp_densplot",
                            "output_duprate_exp_denscurve_mqc": "dupradar_output_duprate_exp_denscurve_mqc",
                            "output_expression_histogram": "dupradar_output_expression_histogram",
                            "output_intercept_slope": "dupradar_output_intercept_slope"
                        ],
                    toState: { id, output, state ->
                        [ 
                            "dupradar_output_dupmatrix": output.output_dupmatrix,
                            "dupradar_output_dup_intercept_mqc": output.output_dup_intercept_mqc,
                            "dupradar_output_duprate_exp_boxplot": output.output_duprate_exp_boxplot,
                            "dupradar_output_duprate_exp_densplot": output.output_duprate_exp_densplot,
                            "dupradar_output_duprate_exp_denscurve_mqc": output.output_duprate_exp_denscurve_mqc,
                            "dupradar_output_expression_histogram": output.output_expression_histogram,
                            "dupradar_output_intercept_slope": output.output_intercept_slope
                    ]
                        }
                )

                & 

                qualimap.run(
                    fromState: [
                            "input": "bam_input",
                            "gtf": "gtf_input",
                            "output_pdf": "qualimap_output_pdf",
                            "output_format": "output_format",
                            "pr_bases": "pr_bases",
                            "tr_bias": "tr_bias",
                            "algorithm": "algorithm",
                            "sequencing_protocol": "sequencing_protocol",
                            "sorted": "sorted",
                            "java_memory_size": "java_memory_size"
                        ],
                    toState: { id, output, state ->
                        [
                            "qualimap_output_pdf": output.output_pdf
                        ]
                    }
                )
                
            )
        | mix()
        | toSortedList({a, b -> a[0] <=> b[0]})
        | map { list -> 
            def ids = list.collect{it[0]}.unique()
            assert ids.size() == 1
            def id = ids[0]
            def state = list.inject([:]){currentState, tuple -> return currentState + tuple[1] }
            [id, state]
        }

        | niceView()

        | setState(
            ["bamstat_output": "bamstat_output",
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
            "output_duplication_rate_plot": "output_duplication_rate_plot",
            "output_duplication_rate_mapping": "output_duplication_rate_mapping",
            "output_duplication_rate_sequence": "output_duplication_rate_sequence",
            "tin_output_summary": "tin_output_summary",
            "tin_output_metrics": "tin_output_metrics",
            "dupradar_output_dupmatrix": "dupradar_output_dupmatrix",
            "dupradar_output_dup_intercept_mqc": "dupradar_output_dup_intercept_mqc",
            "dupradar_output_duprate_exp_boxplot": "dupradar_output_duprate_exp_boxplot",
            "dupradar_output_duprate_exp_densplot": "dupradar_output_duprate_exp_densplot",
            "dupradar_output_duprate_exp_denscurve_mqc": "dupradar_output_duprate_exp_denscurve_mqc",
            "dupradar_output_expression_histogram": "dupradar_output_expression_histogram",
            "dupradar_output_intercept_slope": "dupradar_output_intercept_slope"]
            
            )

    emit:
        output_ch

}
