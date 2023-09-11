workflow run_wf {

    take: 
        input_ch

    main: 
        output_ch = input_ch 
        | view {"Input: $it"}
        | gunzip.run(auto: [publish: true], fromState: ["input": "fasta"], toState: ["fasta_uncompressed": "output"], args: [output: "genome.fasta"])
        | view {"State: $it"}
        | gunzip.run(auto: [publish: true], fromState: ["input": "gtf"], toState: ["gtf_uncompressed": "output"], args: [output: "genes.gtf"])
        | view {"State: $it"}
        | gunzip.run(auto: [publish: true], fromState: ["input": "gff"], toState: ["gff_uncompressed": "output"], args: [output: "genes.gff"])
        | view {"State: $it"}
        | gffread.run(auto: [publish: true], fromState: ["input": "gff_uncompressed"], toState: ["gtf_uncompressed": "output"], args: [output: "genes.gtf"])
        | gunzip.run(auto: [publish: true], fromState: ["input": "additional_fasta"], toState: ["additional_fasta_uncompressed": "output"], args: [output: "additional_fasta.fasta"])
        | view {"State: $it"}
        | cat_additional_fasta.run(auto: [publish: true], fromState: ["fasta": "fasta_uncompressed", "gtf": "gtf_uncompressed", "additional_fasta": "additional_fasta_uncompressed", "biotype": "biotype"], toState: ["concatenated_fasta": "fasta_output", "concatenated_gtf": "gtf_output"], args: [fasta_output: "concatenated.fasta", gtf_output: "concatenated.gtf"])
        | view {"Output: $it"}
    emit: 
        output_ch
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