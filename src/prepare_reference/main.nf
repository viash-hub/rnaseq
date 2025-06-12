workflow run {
  take: input_channel
  main:
    output_channel = input_channel

      | bgzip.run(
        runIf: { id, state ->
          // TODO: check if the input_genome_fasta is gzipped
          def is_fasta_gzipped = true
          is_fasta_gzipped
        },
        fromState: [
          input: "input_genome_fasta"
        ],
        toState: [
          input_genome_fasta: "output"
        ]
      )

      | setState(
        [
          output_genome_fasta: "input_genome_fasta",
          output_transcriptome_gtf: "input_transcriptome_gtf",
        ]
      )

  emit: output_channel
}