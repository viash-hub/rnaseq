workflow run_wf {
    take: 
        input_ch

    main: 
        output_ch = input_ch

        | toSortedList

        | map { list -> 
            salmon_quant_merged = list.collect{id, state -> state.salmon_quant_results}
            [
                "merged", 
                [ salmon_quant_merged: salmon_quant_merged, 
                gtf: list[1][-1].gtf, 
                gtf_extra_attributes: list[1][-1].gtf_extra_attributes, 
                gtf_group_features: list[1][-1].gtf_group_features ] 
            ]
        }
    
        | salmon_tx2gene.run (
            fromState: [ 
                "salmon_quant_results": "salmon_quant_merged", 
                "gtf_extra_attributes": "gtf_extra_attributes", 
                "gtf": "gtf", 
                "gtf_group_features": "gtf_group_features"
            ],
            toState: [ "tx2gene_tsv": "tsv" ]
        )

        | salmon_tximport.run (
            fromState: [ 
                "salmon_quant_results": "salmon_quant_merged", 
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