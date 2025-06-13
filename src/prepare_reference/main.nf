include { isGzipped } from meta.resources_dir + "/is_gzipped.nf"

workflow run_wf {
  take: input_ch
  main:
    output_ch = input_ch

      | bgzip.run(
        runIf: { id, state ->
            isGzipped(state.input_genome_fasta)
        },
        key: "bgzip_genome_fasta",
        fromState: [ input: "input_genome_fasta" ],
        args: [ decompress: true ],
        toState: [ input_genome_fasta: "output" ]
      )

      | bgzip.run(
        runIf: { id, state ->
            isGzipped(state.input_transcriptome_gtf)
        },
        key: "bgzip_transcriptome_gtf",
        fromState: [ input: "input_transcriptome_gtf" ],
        args: [ decompress: true ],
        toState: [ input_transcriptome_gtf: "output" ]
      )

      | star_genome_generate.run(
        fromState: [
          genome_fasta_files: "input_genome_fasta",
          sjdb_gtf_file: "input_transcriptome_gtf"
        ],
        toState: [ output_star_index: "index" ],
      )

      | salmon_index.run(
        fromState: [
          transcripts: "input_genome_fasta"
        ],
        toState: [ output_salmon_index: "index" ]
      )

      | setState(
        [
          output_genome_fasta: "input_genome_fasta",
          output_transcriptome_gtf: "input_transcriptome_gtf",
          output_star_index: "output_star_index",
          output_salmon_index: "output_salmon_index"
        ]
      )

  emit: output_ch
}