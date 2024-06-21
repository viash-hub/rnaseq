workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch
    | map { id, state ->
      def input = state.fastq_2 ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
      def paired = input.size() == 2
      [ id, state + [ paired: paired, input: input ] ]
    }


    // Infer lib-type for salmon quant
    | map { id, state -> 
      def lib_type = (state.paired) ? 
        (
          (state.strandedness == "forward") ? 
            "ISF" : 
            (
              (state.strandedness == "reverse") ? "ISR" : "IU"
            )
        ) 
        : (
          (state.strandedness == "forward") ? 
            "SF" : 
            (
              (state.strandedness == "reverse") ? "SR" : "U"
            )
        ) 
      [ id, state + [lib_type: lib_type] ]
    }
    
    // Count reads from BAM alignments using Salmon
    | salmon_quant.run ( 
        runIf: { id, state -> state.pseudo_aligner == 'salmon' },
        fromState: { id, state ->
          def unmated_reads = !state.paired ? state.subsampled_fastq_1 : null
          def mates1 = state.paired ? state.subsampled_fastq_1 : null
          def mates2 = state.paired ? state.subsampled_fastq_2 : null
          [ unmated_reads: unmated_reads,
            mates1: state.fastq1, 
            mates2: state.fastq2, 
            targets: state.transcript_fasta, 
            gene_map: state.gtf, 
            index: state.pseudo_index,
            lib_type: state.lib_type ]
        },
        toState: [ 
          "quant_results_dir": "output",
          "quant_results_file": "quant_results" 
        ]
    )

    | kallisto_quant.run ( 
        runIf: { id, state -> state.pseudo_aligner == 'kallisto'},
        fromState: [
          "input": "input", 
          "gtf": "gtf", 
          "index": "pseudo_index",
          "fragment_length": "kallisto_quant_fragment_length", 
          "fragment_length_sd": "kallisto_quant_fragment_length_sd"
        ],
        toState: [
          "quant_out_dir": "output", 
          "pseudo_multiqc": "log"
        ]
    )

    | setState (
      [ "pseudo_multiqc": "quant_results", 
        "quant_out_dir": "quant_out_dir",
        "quant_results_file": "quant_results_file" ]
    )

  emit:
    output_ch
}
