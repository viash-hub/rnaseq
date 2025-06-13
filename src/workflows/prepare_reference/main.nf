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
        toState: [ fasta: "output" ]
      )

      | bgzip.run(
        runIf: { id, state ->
            isGzipped(state.input_transcriptome_gtf)
        },
        key: "bgzip_transcriptome_gtf",
        fromState: [ input: "input_transcriptome_gtf" ],
        args: [ decompress: true ],
        toState: [ gtf: "output" ]
      )

      | bgzip.run(
        runIf: { id, state ->
            isGzipped(state.additional_fasta)
        },
        key: "bgzip_additional_fasta",
        fromState: [ input: "additional_fasta" ],
        args: [ decompress: true ],
        toState: [ additional_fasta: "output" ]
      )

      | cat_additional_fasta.run(
        runIf { id, state -> state.additional_fasta },
        fromState: [
          fasta: "fasta",
          gtf: "gtf",
          additional_fasta: "additional_fasta",
          biotype: "biotype"
        ],
        toState: [
          fasta: "fasta_output",
          gtf: "gtf_output"
        ]
      )

      | star_genome_generate.run(
        fromState: [
          genome_fasta_files: "fasta",
          sjdb_gtf_file: "gtf"
        ],
        toState: [ output_star_index: "index" ],
      )

      | salmon_index.run(
        fromState: [
          transcripts: "fasta"
        ],
        toState: [ output_salmon_index: "index" ]
      )

      | setState(
        [
          output_genome_fasta: "fasta",
          output_transcriptome_gtf: "gtf",
          output_star_index: "output_star_index",
          output_salmon_index: "output_salmon_index"
        ]
      )

  emit: output_ch
}