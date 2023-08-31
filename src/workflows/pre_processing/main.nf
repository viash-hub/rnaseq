// Note: some helper functionality is added at the end of this file

workflow run_wf {
  take:
    input_ch

  main:
        
    output_ch = input_ch
    | view { viewEvent(it) }
    // Parse input channel and convert to either paired or unpaired input
    | map{ id, state ->
        (existsInDict(state, "fastq_2"))
          ? [ id, state + [ paired: true, input: [ state.fastq_1, state.fastq_2 ] ] ]
          : [ id, state + [ paired: false, input: [ state.fastq_1 ] ] ]
    }
    // [ id, [ input: ... ]]
    | fastqc.run(
      auto: [publish: true], 
      fromState: [ "paired", "input" ],
      toState: ["fastqc_report": "output"]
    )
    // [ id, [ input: ..., fastqc_report: ... ] ]
    | view { viewEvent(it) }
    // | umitools_extract.run(
    //   auto: [publish: true], 
    //   fromState: ["paired", "input", "umitools_bc_pattern", "umitools_bc_pattern2"],
    //   toState: ["umi_extract_output": "output", "umi_extract_r2_output": "r2_output"]
    // )
    // | view { "Input: $it" }
    // | trimgalore.run(
    //   auto: [publish: true], 
    //   fromState: ["umi_extract_output", "umi_extract_r2_output"],
    //   toState: ["trimgalore_output": "output"]
    // )
    // // [ id, [ input: ..., fastqc_report: ..., umi_extract_output: ... ] ]
    // | view { "Output: $it" }

    // println(params.skip_fastqc)
    // println(params.with_umi)
    // println(params.skip_umi_extract)
    // println(params.skip_trimming)

    // fastqc_ch = input_ch
    // if(!params.skip_fastqc) {
    //   fastqc_ch = input_ch
    //   | view { "Input: $it" }
    //   // [ id, [ input: ... ]]
    //   | fastqc.run(
    //     auto: [publish: true], 
    //     fromState: ["input"],
    //     toState: ["fastqc_report": "output"]
    //   )
    //   | view { "Output: $it" }
    // }

    // umi_extract_ch = fastqc_ch
    // if (params.with_umi && !params.skip_umi_extract) {
    //   umi_extract_ch = fastqc_ch 
    //   | view { "Input: $it" }
    //   // [ id, [ input: ..., fastqc_report: ... ] ]
    //   | umitools_extract.run(
    //     auto: [publish: true], 
    //     fromState: ["input", "umitools_bc_pattern"],
    //     toState: ["umi_extract_output": "output"]
    //   )
    //   | view { "Output: $it" }
    // }

    // trimgalore_ch = umi_extract_ch
    // if (!params.skip_trimming) {
    //   trimgalore_ch = umi_extract_ch
    //   | view { "Input: $it" }
    //   | trimgalore.run(
    //     auto: [publish: true], 
    //     fromState: ["umi_extract_output"],
    //     toState: ["trimgalore_output": "output"]
    //   )
    //   // [ id, [ input: ..., fastqc_report: ..., umi_extract_output: ... ] ]
    //   | view { "Output: $it" }
    // }

    // output_ch = trimgalore_ch

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

// Parse a general Channel event to:
// - render proper Paths
// Modified because `rootDir` is not defined here by default
def parseEvent(root) {
  if (root instanceof List) {
    root.collect{ parseEvent(it) }
  } else if (root instanceof Map) {
    root.collectEntries{ [ (it.key): parseEvent(it.value) ] }
  } else {
    root.toString()
  } 
}

import org.yaml.snakeyaml.Yaml
import org.yaml.snakeyaml.DumperOptions

// Function to be used in Nextflow's view operator
def viewEvent(element) {
  dumperOptions = new DumperOptions()
  dumperOptions.setPrettyFlow(true)
  dumperOptions.setDefaultFlowStyle(DumperOptions.FlowStyle.BLOCK)
  Yaml yaml = new Yaml(dumperOptions)
  yaml.dump(parseEvent(element))
}

// Helper function
def existsInDict(dict, key) {
  return dict.containsKey(key) && dict[key] != ""
}

