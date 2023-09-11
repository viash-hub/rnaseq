// Note: some helper functionality is added at the end of this file

workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch
    // | view {"Input: $it"}
    | view { viewEvent(it) }
    // // Parse input channel and convert to either paired or unpaired input
    | map { id, state ->
         (existsInDict(state, "fastq_2"))
           ? [ id, state + [ "paired": true, "input": [ state.fastq_1, state.fastq_2 ] ] ]
           : [ id, state + [ "paired": false, "input": [ state.fastq_1 ] ] ] }
    | view { viewEvent(it) }
    | map { id, state ->
         (state.gencode)
           ? [ id, state + [ "biotype": "gene_type" ] ]
           : [ id, state + [ "biotype": state.featurecounts_group_type ] ] }
    | view { viewEvent(it) }
    // create list "prepareToolIndices"
    // | prepare_genome.run (
    //     fromState: ["fasta", "gtf", 'gff', "additional_fasta", "transcript_fasta", "gene_bed", "splicesites", "bbsplit_fasta_list", "star_index", "rsem_index", "salmon_index", "hisat2_index", "bbsplit_index", "gencode", "biotype", "prepareToolIndices"],
    //     toState: ["fasta": "uncompressed_fasta", "gtf", "fai", "gene_bed", "transcript_fasta", "chrom_sizes", "splicesites". "bbsplit_index", "star_index", "rsem_index", "hisat2_index", "salmon_index" ]
    // )
    // | pre_processing.run()
    // | view { viewEvent(it) }
    // | genome_alignment.run()
    // | view { viewEvent(it) }
    // | pseudo_alignment.run()
    // | view { viewEvent(it) }
    // | post_processing.run()
    // | view { viewEvent(it) }
    // | final_qc.run()
    // | view { viewEvent(it) }
    // | view { "Output: $it" }

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
