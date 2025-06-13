nextflow.enable.dsl=2

include { prepare_reference } from params.rootDir + "/target/nextflow/prepare_reference/main.nf"
params.resources_test = params.rootDir + "/resources_test/minimal_test/reference/"

workflow test_wf {
  resources_test_file = file(params.resources_test)

  output_ch = Channel.fromList(
    [
      [
        id: "test",
        input_genome_fasta: resources_test_file.resolve("genome.fasta"),
        input_transcriptome_gtf: resources_test_file.resolve("genes.gtf.gz"),
        output_genome_fasta: "genome.fasta",
        output_transcriptome_gtf: "transcriptome.gtf",
        output_star_index: "star_index",
        output_salmon_index: "salmon_index"
      ]
    ]
  )

  | map{ state -> [state.id, state] }

  | prepare_reference

  | view { output ->
    assert output.size() == 2: "Output should contain two elements: [id, state]"

    assert output[0] == "test": "Output ID should be same as input ID"

    assert output[1].output_star_index.isDirectory() && output[1].output_star_index.list().size() > 0: 
      "STAR index directory should exist and not be empty"
    assert output[1].output_salmon_index.isDirectory() && output[1].output_salmon_index.list().size() > 0:
      "SALMON index directory should exist and not be empty"

    "Output: $output"
    }

}