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

    // Sub-sample FastQ files and pseudo-align with Salmon to auto-infer strandedness
    | fq_subsample.run (
      runIf: { id, state -> state.strandedness == 'auto' }, 
      fromState: [ 
        "input": "input", 
        "extra_args": "extra_fq_subsample_args" 
      ],
      toState: [
        "fastq_1": "output_1",
        "fastq_2": "output_2"
      ]
    )

    | salmon_quant.run (
      runIf: { id, state -> state.strandedness == 'auto' }, 
      fromState: { id, state ->
        def input = state.paired ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
        [ paired: state.paired, 
          strandedness: state.strandedness, 
          input: state.input, 
          transcript_fasta: state.transcript_fasta, 
          gtf: state.gtf, 
          salmon_index: state.salmon_index ]
      },
        args: [
          "alignment_mode": false, 
          "lib_type": "A", 
          "extra_salmon_quant_args": "--skipQuant"
        ],
        toState: [
          "salmon_quant_output": "output",
          "salmon_json_info": "json_info"
        ]
    )

    | infer_strandedness.run (
      runIf: { id, state -> state.strandedness == 'auto' }, 
      fromState: [ "json": "salmon_json_info" ], 
      toState: [ "strandedness": "strandedness"], 
    )

    // perform QC on input fastq files
    | fastqc.run (
      runIf: { id, state -> !state.skip_qc && !state.skip_fastqc },
      fromState: [ 
        "paired": "paired",
        "input": "input" 
      ],
      toState: ["fastqc_report": "output"]
    )

    // extract UMIs from fastq files and discard read 1 or read 2 if required
    | umitools_extract.run (
      runIf: { id, state -> state.with_umi && !state.skip_umi_extract },
      fromState: { id, state ->
        def input = state.paired ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
        [ paired: state.paired,
        input: state.input, 
        bc_pattern: state.bc_pattern, 
        umi_discard_read: state.umi_discard_read ]
      },
      toState: { id, state -> 
        def paired = (state.umi_discard_read == 1 || state.umi_discard_read == 2) ? false : true
        [ "fastq_1": "fastq_1", 
        "fastq_2": "fastq_2" ]
      }
    )
    
    // trim reads
    | trimgalore.run (
      runIf: { id, state -> !state.skip_trimming },
      fromState: { id, state ->
        def input = state.paired ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
        [ paired: state.paired,
        input: state.input, 
        extra_trimgalore_args: state.extra_trimgalore_args, 
        min_trimmed_reads: state.min_trimmed_reads ]
      },
      toState: [
        "fastq_1": "fastq_1", 
        "fastq_2": "fastq_2",
        "trim_log": "trim_log", 
        "trim_html": "trim_html", 
        "trim_zip": "trim_zip",
      ]
    )

    // TODO: Filter FATQ files based on minimum trimmed read count after adapter trimming
    // | filter { id, state -> state.skip_trimming || state.trim_read_count >= state.min_trimmed_reads }
    
    // TODO: Get list of samples that failed trimming threshold for MultiQC report

    // filter out rRNA reads
    | bbmap_bbsplit.run (
      runIf: { id, state -> !state.skip_bbsplit },
      fromState: { id, state ->
        def input = state.paired ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
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
    | sortmerna.run (
      runIf: { id, state -> state.remove_ribo_rna },
      fromState: { id, state ->
        def input = state.paired ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
        [ paired: state.paired,
          input: state.input,
          ribo_database_manifest: state.ribo_database_manifest ] 
      },
      toState: [
        "fastq_1": "fastq_1", 
        "fastq_2": "fastq_2",
        "sortmerna_log": "sortmerna_log" ] 
    )

    // Clean up state such that the state only contains arguments
    // with `direction: output` in the viash config
    | setState ( 
        [ "fastqc_report": "fastqc_report", 
        "qc_output1": "fastq_1",
        "qc_output2": "fastq_2", 
        "trim_log": "trim_log", 
        "trim_zip": "trim_zip",
        "sortmerna_log": "sortmerna_log" ] )

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