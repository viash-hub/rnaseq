workflow run_wf {
  take:
    input_ch
  main:
    output_ch = input_ch

      | concat_r1.run(
        fromState: [
          input: "input_r1"
        ],
        toState: [
          processed_r1: "output"
        ]
      )

      | concat_r2.run(
        fromState: [
          input: "input_r2"
        ],
        toState: [
          processed_r2: "output"
        ]
      )

      // TODO: add fq linter

      // run fastqc on raw reads
      | fastqc_raw.run(
        fromState: { id, state ->
          [
            input: [state.processed_r1, state.processed_r2]
          ]
        },
        toState: {
          fastqc_raw_zip: "zip"
        }
      )

      // TODO: add fq trimmer (trimgalore or fastp)

      // TODO: run fastqc on trimmed reads
      // TODO: lint again?

      // TODO: remove genome contaminant reads (bbmap_bbsplit)
      // TODO: lint again?

      // TODO: remove ribosomal RNA reads (sortmerna)
      // TODO: lint again?
      
      // TODO: infer strandedness (if not provided)
      | map { id, state ->
        def newState = state + ["strandedness": "forward"]
        [id, newState]
      }

      | setState(
        [
          output_r1: "processed_r1",
          output_r2: "processed_r2",
          output_strandedness: "strandedness"
        ]
      )
  
  emit:
    output_ch
}
