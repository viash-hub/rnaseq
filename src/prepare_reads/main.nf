workflow run_wf {
  take:
    input_ch
  main:
    output_ch = input_ch

      | concat_text.run(
        fromState: [
          input: "input_r1"
        ],
        toState: [
          processed_r1: "output"
        ]
      )

      | concat_text.run(
        fromState: [
          input: "input_r2"
        ],
        toState: [
          processed_r2: "output"
        ]
      )
  
  emit:
    output_ch
}
