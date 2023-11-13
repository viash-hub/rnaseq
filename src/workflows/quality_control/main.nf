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
                            "map_qual": "map_qual"
                        ],
                    toState: { id, output, state ->
                        [
                            "bamstat_output": output.output
                        ]
                    }
                )

                & 
                
                rseqc_inferexperiment.run(
                    fromState: [
                            "input": "bam_input",
                            "refgene": "refgene",
                            "sample_size": "sample_size",
                            "map_qual": "map_qual"
                        ],
                    toState: { id, output, state ->
                        [
                            "strandedness_output": output.output
                        ]
                    }
                )
            )
        
        | mix

        | view

        | setState(
            "bamstat_output": "bamstat_output",
            "strandedness_output": "strandedness_output"

        )

    emit:
        output_ch

}
