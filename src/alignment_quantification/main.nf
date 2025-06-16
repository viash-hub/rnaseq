workflow run_wf {
  take: input_channel
  main:

    // Align reads to genome and transcriptome using STAR
    star_ch = input_channel
      | star_align_reads.run(
        fromState: [
          genome_dir: "input_star_genome_dir",
          input: "input_r1",
          input_r2: "input_r2",
          sjdb_gtf_file: "input_gtf"
        ],
        args: [
          quant_mode: "TranscriptomeSAM",
          twopass_mode: "Basic",
          out_sam_type: ["BAM", "Unsorted"],
          // read_files_command: "zcat", // See https://github.com/viash-hub/biobox/issues/178
          run_rng_seed: 0,
          out_filter_multimap_nmax: 20,
          // --alignSJDBoverhangMin 1 // Argument not exposed by component
          out_sam_attributes: ["NH", "HI", "AS", "NM", "MD"],
          out_sam_strand_field: "intronMotif",
          quant_transcriptome_sam_output: "BanSingleEnd"
        ],
        toState: [
          star_bam_genome: "aligned_reads",
          star_bam_transcriptome: "reads_aligned_to_transcriptome",
          star_junctions: "splice_junctions",
          star_log: "log"
        ]
      )

    // Quantify expression using salmon in alignment-based mode
    salmon_ch = star_ch
      | salmon_quant.run(
        fromState: [
          alignments: "star_bam_transcriptome",
          targets: "input_transcript_fasta",
          gene_map: "input_gtf"
        ],
        args: [lib_type: "A"],
        toState: [
          salmon_quant: "output"
        ]
      )

    // Set output files
    output_channel = star_ch.join(salmon_ch)
    | map { id, star_state, salmon_state ->
      def output_state = [
        output_star_bam_genome: star_state.star_bam_genome,
        output_star_bam_transcriptome: star_state.star_bam_transcriptome,
        output_star_junctions: star_state.star_junctions,
        output_star_log: star_state.star_log,
        output_salmon: salmon_state.salmon_quant
      ]

      [id, output_state]
    }

  emit: output_channel
}
