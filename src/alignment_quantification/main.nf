workflow run_wf {
  take: input_channel
  main:
    output_channel = input_channel

      | star_align_reads.run(
        fromState: [
          genome_dir: "input_star_genome_dir",
          input: "input_r1",
          input_r2: "input_r2"
        ],
        toState: [
          output_star_bam: "aligned_reads",
          output_star_junctions: "splice_junctions",
          output_star_log: "log"
        ]
      )

      | setState(
        [
          output_star_bam: "output_star_bam",
          output_star_junctions: "output_star_junctions",
          output_star_log: "output_star_log"
        ]
      )

  emit: output_channel
}
