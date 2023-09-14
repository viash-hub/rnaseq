workflow run_wf {

    take: 
        input_ch

    main: 
        output_ch = input_ch 
        | view {"Input: $it"}

        // decompress fasta
        | gunzip.run(auto: [publish: true], fromState: ["input": "fasta"], toState: ["fasta_uncompressed": "output"], key: "gunzip_fasta")//, args: [output: "genome.fasta"])
        | view {"State: $it"}

        // decompress gtf
        | gunzip.run(auto: [publish: true], fromState: ["input": "gtf"], toState: ["gtf_uncompressed": "output"], key: "gunzip_gtf")
        | view {"State: $it"}

        // decompress gff
        | gunzip.run(auto: [publish: true], fromState: ["input": "gff"], toState: ["gff_uncompressed": "output"], key: "gunzip_gff")
        | view {"State: $it"}

        // gff to gtf
        // | gffread.run(auto: [publish: true], fromState: ["input": "gff_uncompressed"], toState: ["gtf_uncompressed": "output"], args: [output: "genes.gtf"])
        // | view {"State: $it"}

        // decompress additional fasta
        | gunzip.run(auto: [publish: true], fromState: ["input": "additional_fasta"], toState: ["additional_fasta_uncompressed": "output"], key: "gunzip_additional_fasta")
        | view {"State: $it"}

        // concatenate additional fasta
        | cat_additional_fasta.run(auto: [publish: true], fromState: ["fasta": "fasta_uncompressed", "gtf": "gtf_uncompressed", "additional_fasta": "additional_fasta_uncompressed", "biotype": "biotype"], toState: ["concatenated_fasta": "fasta_output", "concatenated_gtf": "gtf_output"], key: "cat_additional") 

        // decompress bed file
        // | view {"State: $it"}
        // | gunzip.run(auto: [publish: true], fromState: ["input": "gene_bed"], toState: ["gene_bed_uncompressed": "output"], args: [output: "genes.bed"])
        | view {"State: $it"}

        // gtf to bed -> not working?
        // | gtf2bed.run(auto: [publish: true], fromState: ["gtf": "concatenated_gtf"], toState: ["gene_bed_uncompressed": "bed_output"]) //, args: [bed_output: "genes.bed"])
        // | view {"State: $it"}

        // decompress transcript fasta
        | gunzip.run(auto: [publish: true], fromState: ["input": "transcript_fasta"], toState: ["transcript_fasta_uncompressed": "output"], key: "transcript_fasta")
        | view {"State: $it"}

        // preprocess transcripts fasta if gtf is in gencode format
        | preprocess_transcripts_fasta.run(auto: [publish: true], fromState: ["transcript_fasta": "transcript_fasta_uncompressed"], toState: ["transcript_fasta_fixed": "fixed_fasta"], key: "transcript_fixed")
        | view {"State: $it"}

        // filter gtf for genes in genome
        | gtf_gene_filter.run(auto: [publish: true], fromState: ["fasta": "concatenated_fasta", "gtf": "concatenated_gtf"], toState: ["filtered_gtf": "filtered_gtf"])
        | view {"State: $it"}

        // prepare reference for RSEM
        // | rsem_prepare_reference.run(auto: [publish: true], fromState: ["fasta": "concatenated_fasta", "gtf": "filtered_gtf"], toState: ["transcripts_fasta": "transcript_fasta"], key: "rsem_ref")
        // | view {"State: $it"}

        // chromosome size and fai index
        | getchromsize.run(auto: [publish: true], fromState: ["fasta": "concatenated_fasta", "gtf": "filtered_gtf"], toState: ["sizes": "sizes", "fai": "fai"], key: "chromsizes")
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