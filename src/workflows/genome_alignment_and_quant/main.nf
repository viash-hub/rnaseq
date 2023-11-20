// Note: some helper functionality is added at the end of this file

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
          "star_ignore_sjdbgtf" 
        ],
        toState: [
          "star_alignment": "output", 
          "genome_bam": "star_align_bam", 
          "transcriptome_bam": "star_align_bam_transcriptome", 
          "star_multiqc": "log_final" 
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
          "bam_csi_index": "bam_csi_index" 
        ],
        toState: [ 
          "genome_bam_bai": "output_bai", 
          "genome_bam_csi": "output_csi"
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
          "fasta": "fasta" 
        ],
        toState: ["genome_bam_stats": "output"],
        key: "genome_stats"
    )
    | samtools_flagstat.run (
        fromState: [
          "bam": "genome_bam_sorted", 
          "bai": "genome_bam_index", 
          "fasta": "fasta" 
        ],
        toState: ["genome_bam_flagstat": "output"],
        key: "genome_flagstat"
    )
    | samtools_idxstats.run(
        fromState: [
          "bam": "genome_bam_sorted", 
          "bai": "genome_bam_index", 
          "fasta": "fasta" 
        ],
        toState: ["genome_bam_idxstats": "output"],
        key: "genome_idxstats"
    )

    // TRANSCRIPTOME BAM
    | samtools_sort.run (
        fromState: ["input": "transcriptome_bam"],
        toState: ["transcriptome_bam_sorted": "output"],
        key: "transcriptome_sorted"
    )
    | samtools_index.run (
        fromState: [
          "bam": "transcriptome_bam_sorted", 
          "bam_csi_index": "bam_csi_index" 
        ],
        toState: [
          "transcriptome_bam_bai": "bai", 
          "transcriptome_bam_csi": "csi" 
        ],
        key: "transcriptome_sorted", 
    )
    | map { id, state -> 
      def transcriptome_bam_index = state.transcriptome_bam_bai ? state.transcriptome_bam_bai : state.transcriptome_bam_csi
      [ id, state + [ transcriptome_bam_index: transcriptome_bam_index ] ]
    }
    | samtools_stats.run (
        fromState: [
          "bam": "transcriptome_bam_sorted", 
          "bai": "transcriptome_bam_index" 
        ],
        toState: ["transcriptome_bam_stats": "output"],
        key: "transcriptome_stats"
    )
    | samtools_flagstat.run (
        fromState: [
          "bam": "transcriptome_bam_sorted", 
          "bai": "transcriptome_bam_index" 
        ],
        toState: ["transcriptome_bam_flagstat": "output"],
        key: "transcriptome_flagstat"
    )
    | samtools_idxstats.run(
        fromState: [
          "bam": "transcriptome_bam_sorted", 
          "bai": "transcriptome_bam_index" 
        ],
        toState: ["transcriptome_bam_idxstats": "output"],
        key: "transcriptome_idxstats"
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
          "get_output_stats": "umi_dedup_stats" 
        ],
        toState: ["genome_bam_sorted": "output_bam"],
        key: "genome_deduped"
    )
    | samtools_index.run (
        runIf: { id, state -> state.with_umi },
        fromState: [
          "bam": "genome_bam_sorted", 
          "bam_csi_index": "bam_csi_index"
        ],
        toState: [ 
          "genome_bam_bai": "output_bai", 
          "genome_bam_csi": "output_csi"
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
          "fasta": "fasta" 
        ],
        toState: ["genome_bam_stats": "output"],
        key: "genome_deduped_stats"
    )
    | samtools_flagstat.run (
        runIf: { id, state -> state.with_umi },
        fromState: [
          "bam": "genome_bam_sorted", 
          "bai": "genome_bam_index", 
          "fasta": "fasta" 
        ],
        toState: ["genome_bam_flagstat": "output"],
        key: "genome_deduped_flagstat"
    )
    | samtools_idxstats.run(
        runIf: { id, state -> state.with_umi },
        fromState: [
          "bam": "genome_bam_sorted", 
          "bai": "genome_bam_index", 
          "fasta": "fasta" 
        ],
        toState: ["genome_bam_idxstats": "output"],
        key: "genome_deduped_idxstats"
    )

    // Deduplicate transcriptome BAM file
    | umitools_dedup.run ( 
        runIf: { id, state -> state.with_umi },
        fromState: [
          "paired": "paired", 
          "bam": "transcriptome_bam_sorted", 
          "bai": "transcriptome_bam_index",
          "get_output_stats": "umi_dedup_stats" 
        ],
        toState: ["transcriptome_bam_deduped": "output_bam"],
        key: "transcriptome_deduped"
    )
    | samtools_sort.run (
        runIf: { id, state -> state.with_umi }, 
        fromState: ["input": "transcriptome_bam_deduped"],
        toState: ["transcriptome_bam_sorted": "output"],
        key: "transcriptome_deduped_sorted"
    )
    | samtools_index.run (
      runIf: { id, state -> state.with_umi },
      fromState: [
        "bam": "transcriptome_bam_sorted", 
        "bam_csi_index": "bam_csi_index" 
      ],
      toState: [
        "transcriptome_bam_bai": "bai", 
        "transcriptome_bam_csi": "csi" 
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
          "bam": "transcriptome_bam_sorted", 
          "bai": "transcriptome_bam_index" 
        ],
        toState: ["transcriptome_bam_stats": "output"],
        key: "transcriptome_deduped_stats"
    )
    | samtools_flagstat.run (
        runIf: { id, state -> state.with_umi }, 
        fromState: [
          "bam": "transcriptome_bam_sorted", 
          "bai": "transcriptome_bam_index" 
        ],
        toState: ["transcriptome_bam_flagstat": "output"],
        key: "transcriptome_deduped_flagstat"
    )
    | samtools_idxstats.run(
        runIf: { id, state -> state.with_umi }, 
        fromState: [
          "bam": "transcriptome_bam_sorted", 
          "bai": "transcriptome_bam_index" 
        ],
        toState: ["transcriptome_bam_idxstats": "output"],
        key: "transcriptome_deduped_idxstats"
    ) 

    // Fix paired-end reads in name sorted BAM file
    | umitools_prepareforquant.run (
        runIf: {id, state -> state.with_umi && state.paired},
        fromState: ["bam": "transcriptome_bam_sorted"],
        toState: ["transcriptome_bam_sorted": "output"]
    )

    // Count reads from BAM alignments using Salmon
    | salmon_quant.run ( 
        fromState: [
          "paired": "paired", 
          "strandedness": "strandedness", 
          "input": "transcriptome_bam_sorted", 
          "transcript_fasta": "transcript_fasta", 
          "gtf": "gtf", 
        ],
        args: ["alignment_mode": true],
        toState: ["salmon_quant_results": "output"]
    )

    | setState (
      [ "star_alignment": "star_alignment", 
        "star_multiqc": "star_multiqc", 
        "genome_bam_sorted": "genome_bam_sorted", 
        "genome_bam_index": "genome_bam_index",  
        "genome_bam_stats": "genome_bam_stats", 
        "genome_bam_flagstat": "genome_bam_flagstat", 
        "genome_bam_idxstats": "genome_bam_idxstats", 
        "transcriptome_bam_sorted": "genome_bam_sorted", 
        "transcriptome_bam_index": "transcriptome_bam_index", 
        "transcriptome_bam_stats": "transcriptome_bam_stats", 
        "transcriptome_bam_flagstat": "transcriptome_bam_flagstat", 
        "transcriptome_bam_idxstats": "transcriptome_bam_idxstats",
        "salmon_quant_results": "salmon_quant_results" ]
    )

  emit:
    output_ch
}

// ===============================
// === start of test workflows ===
// ===============================

// workflow test_wf {

//   // allow changing the resources_test dir
//   params.resources_test = params.rootDir + "/resources_test"

//   // or when running from s3: params.resources_test = "s3://openpipelines-data/"
//   testParams = [
//     param_list: [
//     ]
//   ]

//   output_ch =
//     channelFromParams(testParams, config)
//       | view { "Input: $it" }
//       | run_wf
//       | view { output ->
//         assert output.size() == 2 : "outputs should contain two elements; [id, file]"
//         // ...
//         "Output: $output"
//       }
//       | toSortedList()
//       | map { output_list ->
//         assert output_list.size() == 1 : "output channel should contain one event"
//         // ...
//       }
  
// }
