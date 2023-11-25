workflow run_wf {
    take: 
        input_ch

    main: 
        output_ch = input_ch
    
        | salmon_tx2gene.run (
            fromState: [ 
                "salmon_quant_results": "salmon_quant_results", 
                "gtf_extra_attributes": "gtf_extra_attributes", 
                "gtf": "gtf", 
                "gtf_group_features": "gtf_group_features"
            ],
            toState: [ "tx2gene_tsv": "tsv" ]
        )

        | salmon_tximport.run (
            fromState: [ 
                "salmon_quant_results": "salmon_quant_results", 
                "tx2gene_tsv": "tx2gene_tsv" 
            ],
            toState:  [
                "tpm_gene": "tpm_gene",
                "counts_gene": "counts_gene",
                "counts_gene_length_scaled": "counts_gene_length_scaled",
                "counts_gene_scaled": "counts_gene_scaled", 
                "tpm_transcript": "tpm_transcript", 
                "counts_transcript": "counts_transcript"
            ]
        )
                
        | salmon_summarizedexperiment.run (
            fromState: [ 
                "tpm_gene": "tpm_gene",
                "counts_gene": "counts_gene",
                "counts_gene_length_scaled": "counts_gene_length_scaled",
                "counts_gene_scaled": "counts_gene_scaled", 
                "tpm_transcript": "tpm_transcript", 
                "counts_transcript": "counts_transcript", 
                "tx2gene_tsv": "tx2gene_tsv" 
            ],
            toState: [ "salmon_merged_summarizedexperiment": "output" ]
        )

        | setState (
            "tpm_gene": "tpm_gene",
            "counts_gene": "counts_gene",
            "counts_gene_length_scaled": "counts_gene_length_scaled",
            "counts_gene_scaled": "counts_gene_scaled", 
            "tpm_transcript": "tpm_transcript", 
            "counts_transcript": "counts_transcript", 
            "salmon_merged_summarizedexperiment": "salmon_merged_summarizedexperiment"
        )

    emit: 
        output_ch
}