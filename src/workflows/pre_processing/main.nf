
workflow run_wf {
  take:
    input_ch

  main: 

    output_ch = input_ch
      | view { "Input: $it" }
      // [ id, [ input: ... ]]
      | fastqc.run(
        auto: [publish: true], 
        fromState: ["input", "input_r2"],
        toState: ["fastqc_report": "output"]
      )

      // [ id, [ input: ..., fastqc_report: ... ] ]
      | umitools_extract.run(
        auto: [publish: true], 
        fromState: ["input", "input_r2", "umitools_bc_pattern", "umitools_bc_pattern2"],
        toState: ["umi_extract_output": "output", "umi_extract_r2_output": "r2_output"]
      )
      | view { "Input: $it" }
      | trimgalore.run(
        auto: [publish: true], 
        fromState: ["umi_extract_output", "umi_extract_r2_output"],
        toState: ["trimgalore_output": "output"]
      )
      // [ id, [ input: ..., fastqc_report: ..., umi_extract_output: ... ] ]
      | view { "Output: $it" }
      

  emit:
    output_ch
}

// ===============================
// === start of test workflows ===
// ===============================

workflow test_wf {

  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    param_list: [
    ]
  ]

  output_ch =
    channelFromParams(testParams, config)
      | view { "Input: $it" }
      | run_wf
      | view { output ->
        assert output.size() == 2 : "outputs should contain two elements; [id, file]"
        // ...
        "Output: $output"
      }
      | toSortedList()
      | map { output_list ->
        assert output_list.size() == 1 : "output channel should contain one event"
        // ...
      }
  
}
