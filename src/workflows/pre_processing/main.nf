workflow run_wf {
  
  take:
    input_ch

  main:
    output_ch = input_ch
    
    | map { id, state ->
      def input = state.fastq_2 ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
      def paired = input.size() == 2
      [ id, state + [paired: paired, input: input] ]
    }
    // perform QC on input fastq files
    //   input format: [ id, [ paired: ..., input: ... ]]
    | fastqc.run (
      fromState: [ 
        "paired": "paired",
        "input": "input" 
      ],
      toState: ["fastqc_report": "output"]
    )

    // extract UMIs from fastq files
    // input format: [ id, [ paired: ..., input: ..., fastqc_report: ... ] ]
    | umitools_extract.run (
      runIf: {id, state -> state.with_umi},
      fromState: { id, state ->
        def input = state.fastq_2 ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
        [ paired: state.paired,
        input: state.input, 
        bc_pattern: state.bc_pattern ]
      },
      toState: [
        "fastq_1": "fastq_1", 
        "fastq_2": "fastq_2"
      ]
    )

    // trim reads
    // input format: [ id, [ paired: ..., input: ..., fastqc_report: ..., output: ... ] ]
    | trimgalore.run (
      fromState: { id, state ->
        def input = state.fastq_2 ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
        [ paired: state.paired,
        input: state.input ]
      },
      toState: [
        "fastq_1": "fastq_1", 
        "fastq_2": "fastq_2" 
      ]
    )

    // filter out rRNA reads
    // input format: [ id, [ paired: ..., input: ..., fastqc_report: ..., output: ... ] ]
    | bbmap_bbsplit.run (
      fromState: { id, state ->
        def input = state.fastq_2 ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
        [ paired: state.paired,
        input: state.input,
        built_bbsplit_index: state.bbsplit_index,
        bbsplit_fasta_list: state.bbsplit_fasta_list ]
      },
      args: ["only_build_index": false], 
      toState: [
        "fastq_1": "fastq_1", 
        "fastq_2": "fastq_2"
      ]
    )

    // sort reads by rRNA and non-rRNA?
    // input format: [ id, [ paired: ..., input: ..., fastqc_report: ..., output: ... ] ]
    | sortmerna.run (
      // example of skip argument
      // runIf: { id, state -> !state.skip_sort}
      fromState: { id, state ->
        def input = state.fastq_2 ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
        [ paired: state.paired,
          input: state.input,
          ribo_database_manifest: state.ribo_database_manifest ] 
      },
      toState: [
        "fastq_1": "fastq_1", 
        "fastq_2": "fastq_2" ] 
    )

    // Clean up state such that the state only contains arguments
    // with `direction: output` in the viash config
    | setState ( 
        [ "fastqc_report": "fastqc_report", 
        "qc_output1": "fastq_1",
        "qc_output2": "fastq_2" ] )

  emit:
    output_ch
}

// ===============================
// === start of test workflows ===
// ===============================

// workflow test_wf {

//   // allow changing the resources_test dir
//   params.resources_test = params.rootDir + "/resources_test"

//   // or when running from s3: params.resources_test = "s3://openpipelines-data/"
//   testParams = [
//     param_list: [
//     ]
//   ]

//   output_ch =
//     channelFromParams(testParams, config)
//       | view { "Input: $it" }
//       | run_wf
//       | view { output ->
//         assert output.size() == 2 : "outputs should contain two elements; [id, file]"
//         // ...
//         "Output: $output"
//       }
//       | toSortedList()
//       | map { output_list ->
//         assert output_list.size() == 1 : "output channel should contain one event"
//         // ...
//       }
  
// }

// // Parse a general Channel event to:
// // - render proper Paths
// // Modified because `rootDir` is not defined here by default
// def parseEvent(root) {
//   if (root instanceof List) {
//     root.collect{ parseEvent(it) }
//   } else if (root instanceof Map) {
//     root.collectEntries{ [ (it.key): parseEvent(it.value) ] }
//   } else {
//     root.toString()
//   } 
// }

// import org.yaml.snakeyaml.Yaml
// import org.yaml.snakeyaml.DumperOptions

// // Function to be used in Nextflow's view operator
// def viewEvent(element) {
//   dumperOptions = new DumperOptions()
//   dumperOptions.setPrettyFlow(true)
//   dumperOptions.setDefaultFlowStyle(DumperOptions.FlowStyle.BLOCK)
//   Yaml yaml = new Yaml(dumperOptions)
//   yaml.dump(parseEvent(element))
// }

// // Helper function
// def existsInDict(dict, key) {
//   return dict.containsKey(key) && dict[key] != ""
// }

// def getInput()