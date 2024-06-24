workflow run_wf {
    take: 
        input_ch

    main: 
        output_ch = input_ch

        | map { id, state -> 
            def quant_results = state.quant_type == 'kallisto' ? state.kallisto_quant_results : state.salmon_quant_results
            [id, state + [quant_results: quant_results]]
        }
        | tx2gene.run (
            fromState: [ 
                "quant_results": "quant_results", 
                "gtf_extra_attributes": "gtf_extra_attributes", 
                "gtf": "gtf", 
                "gtf_group_features": "gtf_group_features", 
                "quant_type": "quant_type"
            ],
            toState: [ "tx2gene_tsv": "tsv" ]
        )

        | tximport.run (
            fromState: [ 
                "quant_results": "quant_results", 
                "tx2gene_tsv": "tx2gene_tsv", 
                "quant_type": "quant_type"
            ],
            toState:  [
                "tpm_gene": "tpm_gene",
                "counts_gene": "counts_gene",
                "counts_gene_length_scaled": "counts_gene_length_scaled",
                "counts_gene_scaled": "counts_gene_scaled", 
                "tpm_transcript": "tpm_transcript", 
                "counts_transcript": "counts_transcript", 
                "length_gene": "length_gene",
                "length_transcript": "length_transcript"
            ]
        )
        
        | summarizedexperiment.run (
            fromState: [ 
                "tpm_gene": "tpm_gene",
                "counts_gene": "counts_gene",
                "counts_gene_length_scaled": "counts_gene_length_scaled",
                "counts_gene_scaled": "counts_gene_scaled", 
                "tpm_transcript": "tpm_transcript", 
                "counts_transcript": "counts_transcript", 
                "tx2gene_tsv": "tx2gene_tsv"
            ],
            toState: [ "quant_merged_summarizedexperiment": "output" ]
        )

        | setState (
            [ "tpm_gene": "tpm_gene",
            "counts_gene": "counts_gene",
            "counts_gene_length_scaled": "counts_gene_length_scaled",
            "counts_gene_scaled": "counts_gene_scaled", 
            "tpm_transcript": "tpm_transcript", 
            "counts_transcript": "counts_transcript", 
            "quant_merged_summarizedexperiment": "quant_merged_summarizedexperiment" ]
        )

    emit: 
        output_ch
}