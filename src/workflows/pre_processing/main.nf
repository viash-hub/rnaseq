// Note: some helper functionality is added at the end of this file

workflow run_wf {
  take:
    input_ch

  main:

    output_ch = input_ch
    | view { viewEvent(it) }
    // Parse input channel and convert to either paired or unpaired input
    // | map { id, state ->
    //     (existsInDict(state, "fastq_2"))
    //       ? [ id, state + [ "paired": true, "input": [ state.fastq_1, state.fastq_2 ] ] ]
    //       : [ id, state + [ "paired": false, "input": [ state.fastq_1 ] ] ]
    // }
    // [ id, [ paired: ..., input: ... ]]
    | fastqc.run(
      auto: [publish: true], 
      fromState: ["paired", "input"],
      toState: ["fastqc_report": "output"]
    )
    // [ id, [ paired: ..., input: ..., fastqc_report: ... ] ]
    | view { viewEvent(it) }
    | map { id, state ->
        (existsInDict(state, "umitools_bc_pattern2"))
          ? [ id, state + [ bc_pattern: [ state.umitools_bc_pattern, state.umitools_bc_pattern2 ] ] ]
          : [ id, state + [ bc_pattern: [ state.umitools_bc_pattern ] ] ]
    }
    | umitools_extract.run(
      auto: [publish: true], 
      fromState: ["paired", "input", "bc_pattern"],
      toState: ["umi_extract_output": "output"]
    )
    // [ id, [ paired: ..., input: ..., fastqc_report: ..., umi_extract_output: ... ] ]
    | view { viewEvent(it) }
    | trimgalore.run(
      auto: [publish: true], 
      fromState: ["paired": "paired", "input": "umi_extract_output"],
      toState: ["trimgalore_output": "output"]
    )
    // [ id, [ paired: ..., input: ..., fastqc_report: ..., umi_extract_output: ..., trimgalore_output: ... ] ]
    | view { "State: $it" }
    | map { id, state -> [ id, state + [ "only_build_index": false ] ] }    
    | bbmap_bbsplit.run(
      auto: [publish: true], 
      fromState: ["paired": "paired", "input": "trimgalore_output", "built_bbsplit_index": "bbsplit_index", "only_build_index": "only_build_index", "bbsplit_fasta_list": "bbsplit_fasta_list"],
      toState: ["bbsplit_filtered_output": "filtered_output"]
    )
    // [ id, [ paired: ..., input: ..., fastqc_report: ..., umi_extract_output: ..., bbsplit_filtered_output: ... ] ]
    | view { viewEvent(it) }
    | sortmerna.run(
      auto: [publish: true], 
      fromState: ["paired": "paired", "input": "bbsplit_filtered_output", "ribo_database_manifest": "ribo_database_manifest"],
      toState: ["sortmerna_output": "output"]
    )
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
