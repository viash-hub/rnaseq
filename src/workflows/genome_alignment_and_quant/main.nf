// Note: some helper functionality is added at the end of this file

workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch
    | map { id, state ->
      def input = state.fastq_2 ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
      def paired = input.size() == 2
      [ id, state + [paired: paired, input: input] ]
    }
    | star_align.run (
        fromState: [
          "paired", "input", 
          "gtf", "star_index", 
          "extra_star_align_args" ],
        toState: [
          "star_alignment": "output", 
          "genome_bam": "star_align_bam", 
          "transcriptome_bam": "star_align_bam_transcriptome" ]
    )

    // GENOME BAM
    | samtools_sort.run (
        fromState: ["input": "genome_bam"],
        toState: ["genome_bam_sorted": "output"],
        key: "genome_bam_sort"
    )
    | samtools_index.run (
        fromState: [ 
          "input": "genome_bam_sorted", 
          "bam_csi_index": "bam_csi_index" ],
        toState: ["genome_bam_indexed": "output"],
        key: "genome_bam_index"
    )
    | samtools_stats.run (
        fromState: ["input": "genome_bam_indexed"],
        toState: ["genome_bam_stats": "output"],
        key: "genome_bam_stats"
    )
    | samtools_flagstat.run (
        fromState: ["input": "genome_bam_indexed"],
        toState: ["genome_bam_flagstat": "output"],
        key: "genome_bam_flagstat"
    )
    | samtools_idxstats.run(
        fromState: ["input": "genome_bam_indexed"],
        toState: ["genome_bam_idxstats": "output"],
        key: "genome_bam_idxstats"
    )

    // TRANSCRIPTOME BAM
    | samtools_sort.run (
        fromState: ["input": "transcriptome_bam"],
        toState: ["transcriptome_bam_sorted": "output"],
        key: "transcriptome_bam_sort"
    )
    | samtools_index.run (
        fromState: [
          "input": "transcriptome_bam_sorted", 
          "bam_csi_index": "bam_csi_index" ],
        toState: ["transcriptome_bam_indexed": "output"],
        key: "transcriptome_bam_index", 
        debug: true
    )
    | samtools_stats.run (
        fromState: ["input": "transcriptome_bam_indexed"],
        toState: ["transcriptome_bam_stats": "output"],
        key: "transcriptome_bam_stats"
    )
    | samtools_flagstat.run (
        fromState: ["input": "transcriptome_bam_indexed"],
        toState: ["transcriptome_bam_flagstat": "output"],
        key: "transcriptome_bam_flagstat"
    )
    | samtools_idxstats.run(
        fromState: ["input": "transcriptome_bam_indexed"],
        toState: ["transcriptome_bam_idxstats": "output"],
        key: "transcriptome_bam_idxstats"
    )
    
    // Remove duplicate reads from BAM file based on UMIs

    // Deduplicate genome BAM file
    | umitools_dedup.run ( 
        runIf: {id, state -> state.with_umi},
        fromState: ["paired": "paired", "id": "id", "input": "genome_bam_indexed", "get_output_stats": "umi_dedup_stats"],
        toState: ["umitools_genome_deduped": "output"],
        key: "genome_dedup", 
        debug: true
    )
    | samtools_index.run (
        runIf: {id, state -> state.with_umi},
        fromState: ["input": "umitools_genome_deduped", "bam_csi_index": "bam_csi_index"],
        toState: ["genome_bam_indexed": "output"],
        key: "genome_deduped_index"
    )
    | samtools_stats.run (
        runIf: {id, state -> state.with_umi},
        fromState: ["input": "genome_bam_indexed"],
        toState: ["genome_bam_stats": "output"],
        key: "genome_deduped_stats"
    )
    | samtools_flagstat.run (
        runIf: {id, state -> state.with_umi},
        fromState: ["input": "genome_bam_indexed"],
        toState: ["genome_bam_flagstat": "output"],
        key: "genome_deduped_flagstat"
    )
    | samtools_idxstats.run(
        runIf: {id, state -> state.with_umi},
        fromState: ["input": "genome_bam_indexed"],
        toState: ["genome_bam_idxstats": "output"],
        key: "genome_deduped_idxstats"
    )
    
    // Deduplicate transcriptome BAM file
    | umitools_dedup.run (
        runIf: {id, state -> state.with_umi},
        fromState: [
          "paired": "paired", 
          "id": "id", 
          "input": "transcriptome_bam_indexed", 
          "get_output_stats": "umi_dedup_stats"],
        toState: ["umitools_transcriptome_deduped": "output"],
        key: "transcriptome_dedup", 
        debug: true
    )
    | samtools_index.run (
        runIf: {id, state -> state.with_umi},
        fromState: [
          "input": "umitools_transcriptome_deduped", 
          "bam_csi_index": "bam_csi_index"],
        toState: ["transcriptome_bam_indexed": "output"],
        key: "transcriptome_deduped_index"
    )
    | samtools_stats.run (
        runIf: {id, state -> state.with_umi},
        fromState: ["input": "transcriptome_bam_indexed"],
        toState: ["transcriptome_bam_stats": "output"],
        key: "transcriptome_deduped_stats"
    )
    | samtools_flagstat.run (
        runIf: {id, state -> state.with_umi},
        fromState: ["input": "transcriptome_bam_indexed"],
        toState: ["transcriptome_bam_flagstat": "output"],
        key: "transcriptome_deduped_flagstat"
    )
    | samtools_idxstats.run(
        runIf: {id, state -> state.with_umi},
        fromState: ["input": "transcriptome_bam_indexed"],
        toState: ["transcriptome_bam_idxstats": "output"],
        key: "transcriptome_deduped_idxstats"
    )

    // Fix paired-end reads in name sorted BAM file
    | umitools_prepareforquant.run (
        runIf: {id, state -> state.with_umi && state.paired},
        fromState: ["bam": "transcriptome_bam_indexed"],
        toState: ["transcriptome_bam_indexed": "output"]
    )

    // Count reads from BAM alignments using Salmon
    | salmon_quant.run ( 
        fromState: [
          "paired": "paired", 
          "strandedness": "strandedness", 
          "input": "transcriptome_bam_indexed", 
          "transcript_fasta": "transcript_fasta", 
          "gtf": "gtf", 
          "star_index": "star_index"],
        args: ["alignment_mode": true],
        toState: ["salmon_quant_output": "output"]
    )
    | salmon_tx2gene.run (
        fromState: [ 
          "salmon_quant_results": "salmon_quant_output", 
          "gtf_extra_attributes": "gtf_extra_attributes", 
          "gtf": "gtf", 
          "gtf_group_features": "gtf_group_features"],
        toState: ["salmon_tx2gene_tsv": "tsv"]
    )
    | salmon_tximport.run (
        fromState: [ 
          "salmon_quant_results": "salmon_quant_output", 
          "tx2gene_tsv": "salmon_tx2gene_tsv" ],
        toState: ["salmon_tximport": "output"]
    )
    | salmon_summarizedexperiment.run (
        fromState: [ 
          "input": "salmon_tximport", 
          "tx2gene": "salmon_tx2gene_tsv" ],
        toState: [ "salmon_summarizedexperiment": "output" ]
          // "merged_gene_rds": "merged_gene_rds", 
          // "merged_gene_rds_length_scaled": "merged_gene_rds_length_scaled", 
          // "merged_gene_rds_scaled": "merged_gene_rds_scaled", 
          // "merged_transcript_rds": "merged_transcript_rds" ]
    )

    | setState (
      [ "star_alignment": "star_alignment", 
        "genome_bam_indexed": "genome_bam_indexed", 
        "genome_bam_stats": "genome_bam_stats", 
        "genome_bam_flagstat": "genome_bam_flagstat", 
        "genome_bam_idxstats": "genome_bam_idxstats", 
        "transcriptome_bam_indexed": "transcriptome_bam_indexed", 
        "transcriptome_bam_stats": "transcriptome_bam_stats", 
        "transcriptome_bam_flagstat": "transcriptome_bam_flagstat", 
        "transcriptome_bam_idxstats": "transcriptome_bam_idxstats",
        "salmon_quant_results": "salmon_quant_output", 
        "salmon_tx2gene": "salmon_tx2gene_tsv", 
        "salmon_tximport": "salmon_tximport", 
        "salmon_summarizedexperiment": "salmon_summarizedexperiment" ]
        // "merged_gene_rds": "merged_gene_rds", 
        // "merged_gene_rds_length_scaled": "merged_gene_rds_length_scaled", 
        // "merged_gene_rds_scaled": "merged_gene_rds_scaled", 
        // "merged_transcript_rds": "merged_transcript_rds" ]
    )

    | view { "Output: $it" }

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
