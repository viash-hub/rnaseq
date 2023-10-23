// Note: some helper functionality is added at the end of this file

workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch
    | view {"Input: $it"}
    | map { id, state ->
      def input = state.fastq_2 ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
      def paired = input.size() == 2
      [ id, state + [paired: paired, input: input] ]
    }
    | star_align.run (
        auto: [publish: true],
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
        auto: [publish: true],
        fromState: ["input": "genome_bam"],
        toState: ["genome_bam_sorted": "output"],
        key: "genome_bam_sort"
    )
    | samtools_index.run (
        auto: [publish: true],
        fromState: [ 
          "input": "genome_bam_sorted", 
          "bam_csi_index": "bam_csi_index" ],
        toState: ["genome_bam_indexed": "output"],
        key: "genome_bam_index"
    )
    | samtools_stats.run (
        auto: [publish: true],
        fromState: ["input": "genome_bam_indexed"],
        toState: ["genome_bam_stats": "output"],
        key: "genome_bam_stats"
    )
    | samtools_flagstat.run (
        auto: [publish: true],
        fromState: ["input": "genome_bam_indexed"],
        toState: ["genome_bam_flagstat": "output"],
        key: "genome_bam_flagstat"
    )
    | samtools_idxstats.run(
        auto: [publish: true],
        fromState: ["input": "genome_bam_indexed"],
        toState: ["genome_bam_stats": "output"],
        key: "genome_bam_idxstats"
    )

    // TRANSCRIPTOME BAM
    | samtools_sort.run (
        auto: [publish: true],
        fromState: ["input": "transcriptome_bam"],
        toState: ["transcriptome_bam_sorted": "output"],
        key: "transcriptome_bam_sort"
    )
    | samtools_index.run (
        auto: [publish: true],
        fromState: [
          "input": "transcriptome_bam_sorted", 
          "bam_csi_index": "bam_csi_index" ],
        toState: ["transcriptome_bam_indexed": "output"],
        key: "transcriptome_bam_index", 
        debug: true
    )
    | samtools_stats.run (
        auto: [publish: true],
        fromState: ["input": "transcriptome_bam_indexed"],
        toState: ["transcriptome_bam_stats": "output"],
        key: "transcriptome_bam_stats"
    )
    | samtools_flagstat.run (
        auto: [publish: true],
        fromState: ["input": "transcriptome_bam_indexed"],
        toState: ["transcriptome_bam_flagstat": "output"],
        key: "transcriptome_bam_flagstat"
    )
    | samtools_idxstats.run(
        auto: [publish: true],
        fromState: ["input": "transcriptome_bam_indexed"],
        toState: ["transcriptome_bam_stats": "output"],
        key: "transcriptome_bam_idxstats"
    )
    
    // Remove duplicate reads from BAM file based on UMIs

    // Deduplicate genome BAM file
    | umitools_dedup.run ( 
        runIf: {id, state -> state.with_umi},
        auto: [publish: true],
        fromState: ["paired": "paired", "id": "id", "input": "genome_bam_indexed", "get_output_stats": "umi_dedup_stats"],
        toState: ["umitools_genome_deduped": "output"],
        key: "genome_dedup", 
        debug: true
    )
    | samtools_index.run (
        runIf: {id, state -> state.with_umi},
        auto: [publish: true],
        fromState: ["input": "umitools_genome_deduped", "bam_csi_index": "bam_csi_index"],
        toState: ["genome_bam_indexed": "output"],
        key: "genome_deduped_index"
    )
    | samtools_stats.run (
        runIf: {id, state -> state.with_umi},
        auto: [publish: true],
        fromState: ["input": "genome_bam_indexed"],
        toState: ["genome_deduped_stats": "output"],
        key: "genome_deduped_stats"
    )
    | samtools_flagstat.run (
        auto: [publish: true],
        runIf: {id, state -> state.with_umi},
        fromState: ["input": "genome_bam_indexed"],
        toState: ["genome_deduped_flagstat": "output"],
        key: "genome_deduped_flagstat"
    )
    | samtools_idxstats.run(
        auto: [publish: true],
        runIf: {id, state -> state.with_umi},
        fromState: ["input": "genome_bam_indexed"],
        toState: ["genome_deduped_stats": "output"],
        key: "genome_deduped_idxstats"
    )
    
    // Deduplicate transcriptome BAM file
    | umitools_dedup.run (
        auto: [publish: true],
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
        auto: [publish: true],
        runIf: {id, state -> state.with_umi},
        fromState: [
          "input": "umitools_transcriptome_deduped", 
          "bam_csi_index": "bam_csi_index"],
        toState: ["transcriptome_bam_indexed": "output"],
        key: "transcriptome_deduped_index"
    )
    | samtools_stats.run (
        auto: [publish: true],
        runIf: {id, state -> state.with_umi},
        fromState: ["input": "transcriptome_bam_indexed"],
        toState: ["transcriptome_deduped_stats": "output"],
        key: "transcriptome_deduped_stats"
    )
    | samtools_flagstat.run (
        auto: [publish: true],
        runIf: {id, state -> state.with_umi},
        fromState: ["input": "transcriptome_bam_indexed"],
        toState: ["transcriptome_deduped_flagstat": "output"],
        key: "transcriptome_deduped_flagstat"
    )
    | samtools_idxstats.run(
        auto: [publish: true],
        runIf: {id, state -> state.with_umi},
        fromState: ["input": "transcriptome_bam_indexed"],
        toState: ["transcriptome_deduped_stats": "output"],
        key: "transcriptome_deduped_idxstats"
    )

    // Fix paired-end reads in name sorted BAM file
    | umitools_prepareforquant.run (
        auto: [publish: true],
        runIf: {id, state -> state.with_umi && state.paired},
        fromState: ["bam": "transcriptome_bam_indexed"],
        toState: ["transcriptome_deduped_indexed": "output"]
    )

    // Count reads from BAM alignments using Salmon
    | salmon_quant.run ( 
        auto: [publish: true],
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
        auto: [publish: true],
        fromState: [ 
          "salmon_quant_results": "salmon_quant_output", 
          "gtf_extra_attributes": "gtf_extra_attributes", 
          "gtf": "gtf", 
          "gtf_group_features": "gtf_group_features"],
        toState: ["salmon_tx2gene_tsv": "tsv"]
    )


    | view { "Output: $it" }

  emit:
    output_ch
}

// ===============================
// === start of test workflows ===
// ===============================

workflow test_wf {

  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    param_list: [
    ]
  ]

  output_ch =
    channelFromParams(testParams, config)
      | view { "Input: $it" }
      | run_wf
      | view { output ->
        assert output.size() == 2 : "outputs should contain two elements; [id, file]"
        // ...
        "Output: $output"
      }
      | toSortedList()
      | map { output_list ->
        assert output_list.size() == 1 : "output channel should contain one event"
        // ...
      }
  
}
