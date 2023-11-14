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
        toState: ["salmon_quant_output": "output"]
    )

    //
    // Filter channels to get samples that passed STAR minimum mapping percentage
    //
    // ch_fail_mapping_multiqc = Channel.empty()
    // if (!params.skip_alignment && params.aligner.contains('star')) {
    //     ch_star_multiqc
    //         .map { meta, align_log -> [ meta ] + WorkflowRnaseq.getStarPercentMapped(params, align_log) }
    //         .set { ch_percent_mapped }

    //     ch_genome_bam
    //         .join(ch_percent_mapped, by: [0])
    //         .map { meta, ofile, mapped, pass -> if (pass) [ meta, ofile ] }
    //         .set { ch_genome_bam }

    //     ch_genome_bam_index
    //         .join(ch_percent_mapped, by: [0])
    //         .map { meta, ofile, mapped, pass -> if (pass) [ meta, ofile ] }
    //         .set { ch_genome_bam_index }

    //     ch_percent_mapped
    //         .branch { meta, mapped, pass ->
    //             pass: pass
    //                 pass_mapped_reads[meta.id] = true
    //                 return [ "$meta.id\t$mapped" ]
    //             fail: !pass
    //                 pass_mapped_reads[meta.id] = false
    //                 return [ "$meta.id\t$mapped" ]
    //         }
    //         .set { ch_pass_fail_mapped }

    //     ch_pass_fail_mapped
    //         .fail
    //         .collect()
    //         .map {
    //             tsv_data ->
    //                 def header = ["Sample", "STAR uniquely mapped reads (%)"]
    //                 WorkflowRnaseq.multiqcTsvFromList(tsv_data, header)
    //         }
    //         .set { ch_fail_mapping_multiqc }
    // }

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
        "salmon_quant_results": "salmon_quant_output" ]
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

//
// Function that parses and returns the alignment rate from the STAR log output
//
public static ArrayList getStarPercentMapped(params, align_log) {
  def percent_aligned = 0
  def pattern = /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/
  align_log.eachLine { line ->
    def matcher = line =~ pattern
    if (matcher) {
        percent_aligned = matcher[0][1].toFloat()
    }
  }
  def pass = false
  if (percent_aligned >= params.min_mapped_reads.toFloat()) {
    pass = true
  }
  return [ percent_aligned, pass ]
}

//
// Create MultiQC tsv custom content from a list of values
//
def multiqcTsvFromList(tsv_data, header) {
  def tsv_string = ""
  if (tsv_data.size() > 0) {
    tsv_string += "${header.join('\t')}\n"
    tsv_string += tsv_data.join('\n')
  }
  return tsv_string
}