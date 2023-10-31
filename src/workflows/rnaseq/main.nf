workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch
    
    // prepare all the necessary files for reference genome
    | prepare_genome.run ( 
        auto: [publish: true], 
        fromState: [
          "fasta": "fasta", 
          "gtf": "gtf", 
          "additional_fasta": "additional_fasta", 
          "transcript_fasta": "transcript_fasta", 
          "bbsplit_fasta_list": "bbsplit_fasta_list", 
          "gencode": "gencode", 
          "biotype": "biotype"],
        toState: [
          "fasta": "uncompressed_fasta", 
          "gtf": "gtf_uncompressed", 
          "transcript_fasta": "transcript_fasta_uncompressed", 
          "fai": "fai", 
          "chrom_sizes": "chrom_sizes", 
          "bbsplit_index": "bbsplit_index_uncompressed", 
          "star_index": "star_index_uncompressed", 
          "salmon_index": "salmon_index_uncompressed", 
          "gene_bed": "gene_bed_uncompressed" ],
        debug: true
    )

    // Parse input channel and convert to either paired or unpaired input
    | map { id, state ->
      def input = state.fastq_2 ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
      def paired = input.size() == 2
      def biotype = state.gencode ? "gene_type" : state.featurecounts_group_type 
      def bc_pattern = state.umitools_bc_pattern2 ? [ state.umitools_bc_pattern, state.umitools_bc_pattern2 ] : [ state.umitools_bc_pattern ] 
      [ id, state + [paired: paired, biotype: biotype, bc_pattern: bc_pattern] ]
    }
    
    // create list "prepareToolIndices"
    
    // Pre-process fastq files
    | pre_processing.run ( 
        auto: [publish: true],
        fromState: [
          "id": "id", 
          "fastq_1": "fastq_1",
          "fastq_2": "fastq_2", 
          "strandedness": "strandedness",
          "with_umi": "with_umi", 
          "bbsplit_index": "bbsplit_index", 
          "bbsplit_fasta_list": "bbsplit_fasta_list", 
          "biotype": "biotype", 
          "bc_pattern": "bc_pattern", 
          "ribo_database_manifest": "ribo_database_manifest"], 
        toState: [ 
          "fastqc_report": "fastqc_report", 
          "fastq_1": "qc_output1",
          "fastq_2": "qc_output2"
        ], 
        debug: true
    )
    | genome_alignment_and_quant.run (
        auto: [publish: true], 
        fromState: [
          "id": "id", 
          "fastq_1": "fastq_1",
          "fastq_2": "fastq_2", 
          "strandedness": "strandedness", 
          "gtf": "gtf", 
          "star_index": "star_index", 
          "extra_star_align_args": "extra_star_align_args", 
          "star_ignore_sjdbgtf": "star_ignore_sjdbgtf",
          "seq_platform": "seq_platform", 
          "seq_center": "seq_center",
          "with_umi": "with_umi", 
          "umi_dedup_stats": "umi_dedup_stats",
          "gtf_group_features": "gtf_group_features",
          "gtf_extra_attributes": "gtf_extra_attributes",
          "salmon_quant_libtype": "salmon_quant_libtype"
        ],
        toState: [
          "star_alignment": "star_alignment", 
          "star_multiqc": "star_multiqc", 
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
          "salmon_summarizedexperiment": "salmon_summarizedexperiment"
        ], 
        debug: true
    )
    // | post_processing_and_qc.run()
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
