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
                "gtf_group_features": "gtf_group_features", 
                "versions": "versions"
            ],
            toState: [ 
                "tx2gene_tsv": "tsv", 
                "versions": "updated_versions" 
            ]
        )

        | salmon_tximport.run (
            fromState: [ 
                "salmon_quant_results": "salmon_quant_results", 
                "tx2gene_tsv": "tx2gene_tsv", 
                "versions": "versions" 
            ],
            toState:  [
                "tpm_gene": "tpm_gene",
                "counts_gene": "counts_gene",
                "counts_gene_length_scaled": "counts_gene_length_scaled",
                "counts_gene_scaled": "counts_gene_scaled", 
                "tpm_transcript": "tpm_transcript", 
                "counts_transcript": "counts_transcript", 
                "versions": "updated_versions"
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
                "tx2gene_tsv": "tx2gene_tsv", 
                "versions": "versions" 
            ],
            toState: [ 
                "salmon_merged_summarizedexperiment": "output", 
                "versions": "updated_versions" 
            ]
        )

        | setState (
            "tpm_gene": "tpm_gene",
            "counts_gene": "counts_gene",
            "counts_gene_length_scaled": "counts_gene_length_scaled",
            "counts_gene_scaled": "counts_gene_scaled", 
            "tpm_transcript": "tpm_transcript", 
            "counts_transcript": "counts_transcript", 
            "salmon_merged_summarizedexperiment": "salmon_merged_summarizedexperiment", 
            "updated_versions": "versions"
        )

    emit: 
        output_ch
}