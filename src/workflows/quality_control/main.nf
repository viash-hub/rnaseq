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
                    toState: {id, output, state -> 
                        [
                            "inner_dist_output_stats": output.output_stats,
                            "inner_dist_output_dist": output.output_dist,
                            "inner_dist_output_freq": output.output_freq,
                            "inner_dist_output_plot": output.output_plot,
                            "inner_dist_output_plot_r": output.output_plot_r
                        ]},
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
