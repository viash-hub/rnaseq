workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch
    | map { id, state ->
      def input = state.fastq_2 ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
      def paired = input.size() == 2
      [ id, state + [ paired: paired, input: input ] ]
    }
    | star_align.run (
        fromState: [
          "input", 
          "gtf", 
          "star_index", 
          "extra_star_align_args", 
          "seq_platform", 
          "seq_center", 
          "star_ignore_sjdbgtf", 
          "versions" 
        ],
        toState: [
          "star_alignment": "output", 
          "genome_bam": "star_align_bam", 
          "transcriptome_bam": "star_align_bam_transcriptome", 
          "star_multiqc": "log_final", 
          "versions": "updated_versions" 
        ]
    )

    // GENOME BAM
    | samtools_sort.run (
        fromState: ["input": "genome_bam"],
        toState: ["genome_bam_sorted": "output"],
        key: "genome_sorted"    )
    | samtools_index.run (
        fromState: [ 
          "input": "genome_bam_sorted", 
          "bam_csi_index": "bam_csi_index", 
          "versions": "versions" 
        ],
        toState: [ 
          "genome_bam_bai": "output_bai", 
          "genome_bam_csi": "output_csi", 
          "versions": "updated_versions"
        ],
        key: "genome_sorted"
    )
    | map { id, state -> 
      def genome_bam_index = state.genome_bam_bai ? state.genome_bam_bai : state.genome_bam_csi
      [ id, state + [ genome_bam_index: genome_bam_index ] ]
    }
    | samtools_stats.run (
        fromState: [
          "bam": "genome_bam_sorted", 
          "bai": "genome_bam_index", 
          "fasta": "fasta", 
          "versions": "versions" 
        ],
        toState: [
          "genome_bam_stats": "output", 
          "versions": "updated_versions"
        ],
        key: "genome_stats"
    )
    | samtools_flagstat.run (
        fromState: [
          "bam": "genome_bam_sorted", 
          "bai": "genome_bam_index", 
          "fasta": "fasta", 
          "versions": "versions" 
        ],
        toState: [
          "genome_bam_flagstat": "output", 
          "versions": "updated_versions"
        ],
        key: "genome_flagstat"
    )
    | samtools_idxstats.run(
        fromState: [
          "bam": "genome_bam_sorted", 
          "bai": "genome_bam_index", 
          "fasta": "fasta", 
          "versions": "versions" 
        ],
        toState: [
          "genome_bam_idxstats": "output", 
          "versions": "updated_versions"
        ],
        key: "genome_idxstats"
    )

    //
    // Remove duplicate reads from BAM file based on UMIs
    // 
    
    // Deduplicate genome BAM file
    | umitools_dedup.run ( 
        runIf: { id, state -> state.with_umi },
        fromState: [
          "paired": "paired", 
          "bam": "genome_bam_sorted", 
          "bai": "genome_bam_index",
          "get_output_stats": "umi_dedup_stats", 
          "versions": "versions" 
        ],
        toState: [
          "genome_bam_sorted": "output_bam", 
          "versions": "updated_versions"
        ],
        key: "genome_deduped"
    )
    | samtools_index.run (
        runIf: { id, state -> state.with_umi },
        fromState: [
          "input": "genome_bam_sorted", 
          "bam_csi_index": "bam_csi_index", 
          "versions": "versions"
        ],
        toState: [ 
          "genome_bam_bai": "output_bai", 
          "genome_bam_csi": "output_csi", 
          "versions": "updated_versions"
        ],
        key: "genome_deduped"
    )
    | map { id, state -> 
      def genome_bam_index = state.genome_bam_bai ? state.genome_bam_bai : state.genome_bam_csi
      [ id, state + [ genome_bam_index: genome_bam_index ] ]
    }
    | samtools_stats.run (
        runIf: { id, state -> state.with_umi },
        fromState: [
          "bam": "genome_bam_sorted", 
          "bai": "genome_bam_index", 
          "fasta": "fasta", 
          "versions": "versions" 
        ],
        toState: [
          "genome_bam_stats": "output", 
          "versions": "updated_versions"
        ],
        key: "genome_deduped_stats"
    )
    | samtools_flagstat.run (
        runIf: { id, state -> state.with_umi },
        fromState: [
          "bam": "genome_bam_sorted", 
          "bai": "genome_bam_index", 
          "fasta": "fasta", 
          "versions": "versions" 
        ],
        toState: [
          "genome_bam_flagstat": "output", 
          "versions": "updated_versions"
        ],
        key: "genome_deduped_flagstat"
    )
    | samtools_idxstats.run(
        runIf: { id, state -> state.with_umi },
        fromState: [
          "bam": "genome_bam_sorted", 
          "bai": "genome_bam_index", 
          "fasta": "fasta", 
          "versions": "versions" 
        ],
        toState: [
          "genome_bam_idxstats": "output", 
          "versions": "updated_versions"
        ],
        key: "genome_deduped_idxstats"
    )

    // Deduplicate transcriptome BAM file

    | samtools_sort.run (
        runIf: { id, state -> state.with_umi },
        fromState: [
          "input": "transcriptome_bam", 
          "versions": "versions"
        ],
        toState: [
          "transcriptome_bam": "output", 
          "versions": "updated_versions"
        ],
        key: "transcriptome_sorted"
    )
    | samtools_index.run (
        runIf: { id, state -> state.with_umi },
        fromState: [
          "input": "transcriptome_bam", 
          "bam_csi_index": "bam_csi_index", 
          "versions": "versions" 
        ],
        toState: [
          "transcriptome_bam_bai": "output_bai", 
          "transcriptome_bam_csi": "output_csi", 
          "versions": "updated_versions" 
        ],
        key: "transcriptome_sorted", 
    )
    | map { id, state -> 
      if (state.with_umi) {
        def transcriptome_bam_index = state.transcriptome_bam_bai ? state.transcriptome_bam_bai : state.transcriptome_bam_csi
        [ id, state + [ transcriptome_bam_index: transcriptome_bam_index ] ]
      } else {
        [ id, state ]
      }
    }
    | samtools_stats.run (
        runIf: { id, state -> state.with_umi },
        fromState: [
          "bam": "transcriptome_bam", 
          "bai": "transcriptome_bam_index", 
          "versions": "versions" 
        ],
        toState: [
          "transcriptome_bam_stats": "output", 
          "versions": "updated_versions"
        ],
        key: "transcriptome_stats"
    )
    | samtools_flagstat.run (
        runIf: { id, state -> state.with_umi },
        fromState: [
          "bam": "transcriptome_bam", 
          "bai": "transcriptome_bam_index", 
          "versions": "versions" 
        ],
        toState: [
          "transcriptome_bam_flagstat": "output", 
          "versions": "updated_versions"
        ],
        key: "transcriptome_flagstat"
    )
    | samtools_idxstats.run(
        runIf: { id, state -> state.with_umi },
        fromState: [
          "bam": "transcriptome_bam", 
          "bai": "transcriptome_bam_index", 
          "versions": "versions" 
        ],
        toState: [
          "transcriptome_bam_idxstats": "output", 
          "versions": "updated_versions"
        ],
        key: "transcriptome_idxstats"
    )    
     
    | umitools_dedup.run ( 
        runIf: { id, state -> state.with_umi },
        fromState: [
          "paired": "paired", 
          "bam": "transcriptome_bam", 
          "bai": "transcriptome_bam_index",
          "get_output_stats": "umi_dedup_stats", 
          "versions": "versions" 
        ],
        toState: [
          "transcriptome_bam_deduped": "output_bam", 
          "versions": "updated_versions"
        ],
        key: "transcriptome_deduped"
    )
    | samtools_sort.run (
        runIf: { id, state -> state.with_umi }, 
        fromState: [
          "input": "transcriptome_bam_deduped", 
          "versions": "versions"
        ],
        toState: [
          "transcriptome_bam": "output", 
          "versions": "updated_versions"
        ],
        key: "transcriptome_deduped_sorted"
    )
    | samtools_index.run (
      runIf: { id, state -> state.with_umi },
      fromState: [
        "input": "transcriptome_bam", 
        "bam_csi_index": "bam_csi_index", 
        "versions": "versions" 
      ],
      toState: [
        "transcriptome_bam_bai": "output_bai", 
        "transcriptome_bam_csi": "output_csi", 
        "versions": "updated_versions" 
      ],
      key: "transcriptome_deduped_sorted", 
    )
    | map { id, state -> 
      def transcriptome_bam_index = state.transcriptome_bam_bai ? state.transcriptome_bam_bai : state.transcriptome_bam_csi
      [ id, state + [ transcriptome_bam_index: transcriptome_bam_index ] ]
    }
    | samtools_stats.run (
        runIf: { id, state -> state.with_umi }, 
        fromState: [
          "bam": "transcriptome_bam", 
          "bai": "transcriptome_bam_index", 
          "versions": "versions" 
        ],
        toState: [
          "transcriptome_bam_stats": "output", 
          "versions": "updated_versions"
        ],
        key: "transcriptome_deduped_stats"
    )
    | samtools_flagstat.run (
        runIf: { id, state -> state.with_umi }, 
        fromState: [
          "bam": "transcriptome_bam", 
          "bai": "transcriptome_bam_index", 
          "versions": "versions" 
        ],
        toState: [
          "transcriptome_bam_flagstat": "output", 
          "versions": "updated_versions"
        ],
        key: "transcriptome_deduped_flagstat"
    )
    | samtools_idxstats.run(
        runIf: { id, state -> state.with_umi }, 
        fromState: [
          "bam": "transcriptome_bam", 
          "bai": "transcriptome_bam_index", 
          "versions": "versions" 
        ],
        toState: [
          "transcriptome_bam_idxstats": "output", 
          "versions": "updated_versions"
        ],
        key: "transcriptome_deduped_idxstats"
    ) 

    // Fix paired-end reads in name sorted BAM file
    | umitools_prepareforquant.run (
        runIf: {id, state -> state.with_umi && state.paired},
        fromState: [
          "bam": "transcriptome_bam", 
          "versions": "versions"
        ],
        toState: [
          "transcriptome_bam": "output", 
          "versions": "updated_versions"
        ]
    )

    // Count reads from BAM alignments using Salmon
    | salmon_quant.run ( 
        fromState: [
          "paired": "paired", 
          "strandedness": "strandedness", 
          "input": "transcriptome_bam", 
          "transcript_fasta": "transcript_fasta", 
          "gtf": "gtf", 
          "versions": "versions" 
        ],
        args: ["alignment_mode": true],
        toState: [
          "salmon_quant_results": "output", 
          "versions": "updated_versions"
        ]
    )

    | setState (
      [ "star_alignment": "star_alignment", 
        "star_multiqc": "star_multiqc", 
        "genome_bam_sorted": "genome_bam_sorted", 
        "genome_bam_index": "genome_bam_index",  
        "genome_bam_stats": "genome_bam_stats", 
        "genome_bam_flagstat": "genome_bam_flagstat", 
        "genome_bam_idxstats": "genome_bam_idxstats", 
        "transcriptome_bam": "transcriptome_bam", 
        "transcriptome_bam_index": "transcriptome_bam_index", 
        "transcriptome_bam_stats": "transcriptome_bam_stats", 
        "transcriptome_bam_flagstat": "transcriptome_bam_flagstat", 
        "transcriptome_bam_idxstats": "transcriptome_bam_idxstats",
        "salmon_quant_results": "salmon_quant_results", 
        "updated_versions": "versions" ]
    )

  emit:
    output_ch
}
