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
            def unmated_reads = !state.paired ? state.fastq_1 : null
            def mates1 = state.paired ? state.fastq_1 : null
            def mates2 = state.paired ? state.fastq_2 : null
            [ unmated_reads: unmated_reads,
              mates1: mates1, 
              mates2: mates2, 
              gene_map: state.gtf, 
              index: state.salmon_index,
              lib_type: state.lib_type ]
          },
          toState: [ 
            "quant_out_dir": "output",
            "salmon_quant_results_file": "quant_results" 
          ]
      )

      | map { id, state -> 
        def mod_state = (state.pseudo_aligner == 'salmon') ? state + [pseudo_multiqc: state.quant_out_dir] : state
        [ id, mod_state ]
      }

      | kallisto_quant.run ( 
          runIf: { id, state -> state.pseudo_aligner == 'kallisto'},
          fromState: [
            "input": "input", 
            "paired": "paired",
            "gtf": "gtf", 
            "index": "kallisto_index",
            "fragment_length": "kallisto_quant_fragment_length", 
            "fragment_length_sd": "kallisto_quant_fragment_length_sd"
          ],
          toState: [
            "quant_out_dir": "output", 
            "kallisto_quant_results_file": "quant_results_file", 
            "pseudo_multiqc": "log"
          ]
      )

      | map { id, state -> 
        def mod_state = state.findAll { key, value -> value instanceof java.nio.file.Path && value.exists() }
        [ id, mod_state ]
      }

      | setState (
        [ "pseudo_multiqc": "quant_results", 
          "quant_out_dir": "quant_out_dir",
          "salmon_quant_results_file": "salmon_quant_results_file",
          "kallisto_quant_results_file": "kallisto_quant_results_file" ]
      )

  emit:
    output_ch
}
