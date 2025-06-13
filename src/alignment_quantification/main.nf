workflow run_wf {
  take: input_channel
  main:
    output_channel = input_channel

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
          // --alignSJDBoverhangMin 1
          out_sam_attributes: ["NH", "HI", "AS", "NM", "MD"],
          out_sam_strand_field: "intronMotif",
          quant_transcriptome_sam_output: "BanSingleEnd"
        ],
        toState: [
          output_star_bam_genome: "aligned_reads",
          output_star_bam_transcriptome: "reads_aligned_to_transcriptome",
          output_star_junctions: "splice_junctions",
          output_star_log: "log"
        ]
      )

      | salmon_quant.run(
        fromState: [
          alignments: "output_star_bam_transcriptome",
          targets: "input_transcript_fasta",
          gene_map: "input_gtf"
        ],
        args: [lib_type: "A"],
        toState: [
          output_salmon: "output"
        ]
      )

      | setState(
        [
          output_star_bam_genome: "output_star_bam_genome",
          output_star_bam_transcriptome: "output_star_bam_transcriptome",
          output_star_junctions: "output_star_junctions",
          output_star_log: "output_star_log",
          output_salmon: "output_salmon"
        ]
      )

  emit: output_channel
}
