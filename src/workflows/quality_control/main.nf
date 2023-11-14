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
                        },
                    auto: [ publish: true ]

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
                        },
                    auto: [ publish: true ]
                )

                & 
                
                rseqc_innerdistance.run(
                    // runIf: {id, state -> state.paired},
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
                        ]},
                    auto: [ publish: true ]
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

                    },
                    auto: [ publish: true ]
                )
            )

        | mix

        | view

        // | setState(
        //     "bamstat_output": "bamstat_output",
        //     "strandedness_output": "strandedness_output",
        // )

    emit:
        output_ch

}
