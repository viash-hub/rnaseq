// Note: some helper functionality is added at the end of this file

workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch
    | view {"Input: $it"}
    | map { id, state ->
         (existsInDict(state, "fastq_2"))
           ? [ id, state + [ "paired": true, "input": [ state.fastq_1, state.fastq_2 ] ] ]
           : [ id, state + [ "paired": false, "input": [ state.fastq_1 ] ] ] 
    }
    | view { viewEvent(it) }
    | star_align.run (
        auto: [publish: true],
        fromState: ["paired", "input", "gtf", "star_index", "extra_star_align_args"],
        toState: ["star_alignment": "output"]
    )
    | view { viewEvent(it) }
    | samtools_sort.run (
        auto: [publish: true],
        fromState: ["input": "star_alignment"],
        toState: ["star_alignment_sorted": "output"]
    )
    | samtools_index.run (
        auto: [publish: true],
        fromState: ["input": "star_alignment_sorted", "bam_csi_index": "bam_csi_index"],
        toState: ["star_alignment_indexed": "output"]
    )
    | samtools_stats.run (
        auto: [publish: true],
        fromState: ["input": "star_alignment_indexed"],
        toState: ["star_alignment_stats": "output"]
    )
    | samtools_flagstat.run (
        auto: [publish: true],
        fromState: ["input": "star_alignment_indexed"],
        toState: ["star_alignment_flagstat": "output"]
    )
    | samtools_idxstats.run(
        auto: [publish: true],
        fromState: ["input": "star_alignment_indexed"],
        toState: ["star_alignment_stats": "output"]
    )
    | view { viewEvent(it) }
    | umitools_dedup.run (
        auto: [publish: true],
        fromState: ["paired": "paired", "input": "star_alignment_indexed"],
        toState: ["umitools_deduped": "output"]
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
