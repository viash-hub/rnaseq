workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch

    | map { id, state ->
      def biotype = state.gencode ? "gene_type" : state.featurecounts_group_type 
      def bc_pattern = state.umitools_bc_pattern2 ? [ state.umitools_bc_pattern, state.umitools_bc_pattern2 ] : [ state.umitools_bc_pattern ] 
      [ id, state + [biotype: biotype, bc_pattern: bc_pattern] ]
    } 

    // prepare all the necessary files for reference genome
    | prepare_genome.run ( 
        auto: [publish: true], 
        fromState: [
          "fasta": "fasta", 
          "gtf": "gtf", 
          "gff": "gff"
          "additional_fasta": "additional_fasta", 
          "transcript_fasta": "transcript_fasta", 
          "gene_bed": "gene_bed",
          "splicesites": "splicesites",
          "bbsplit_fasta_list": "bbsplit_fasta_list", 
          "star_index": "star_index", 
          "rsem_index": "rsem_index",
          "salmon_index": "salmon_index",
          "hisat2_index": "hisat2_index",
          "bbsplit_index": "bbsplit_index",
          "gencode": "gencode", 
          "biotype": "biotype" 
        ],
        toState: [
          "fasta": "uncompressed_fasta", 
          "gtf": "gtf_uncompressed", 
          "transcript_fasta": "transcript_fasta_uncompressed", 
          "fai": "fai", 
          "chrom_sizes": "chrom_sizes", 
          "bbsplit_index": "bbsplit_index_uncompressed", 
          "star_index": "star_index_uncompressed", 
          "salmon_index": "salmon_index_uncompressed", 
          "gene_bed": "gene_bed_uncompressed" 
        ],
        debug: true
    )

    // Check if contigs in genome fasta file > 512 Mbp
    | check_contig_size.run ( 
        fromState: ["input": "fai"],
        debug: true
    )

    // Concatenate FastQ files from same sample if required
    | cat_fastq.run (
        runIf: {id, state -> state.fastq_1.size() > 1 || state.fastq_2.size() > 1}, 
        fromState: [
          "read_1": "fastq_1", 
          "read_2": "fastq_2"
        ], 
        toState: [ 
          "fastq_1": "fastq_1",
          "fastq_2": "fastq_2"
        ], 
        debug: true
    )
    
    // Pre-process fastq files
    | pre_processing.run ( 
        auto: [publish: true],
        fromState: [
          "id": "id", 
          "fastq_1": "fastq_1",
          "fastq_2": "fastq_2", 
          "strandedness": "strandedness",
          "transcript_fasta": "transcript_fasta", 
          "gtf": "gtf",
          "with_umi": "with_umi", 
          "bbsplit_index": "bbsplit_index", 
          "bbsplit_fasta_list": "bbsplit_fasta_list", 
          "bc_pattern": "bc_pattern", 
          "ribo_database_manifest": "ribo_database_manifest" 
        ], 
        toState: [ 
          "fastqc_report": "fastqc_report", 
          "fastq_1": "qc_output1",
          "fastq_2": "qc_output2"
        ], 
        debug: true
    )
    // | genome_alignment_and_quant.run (
    //     auto: [publish: true], 
    //     fromState: [
    //       "id": "id", 
    //       "fastq_1": "fastq_1",
    //       "fastq_2": "fastq_2", 
    //       "strandedness": "strandedness", 
    //       "gtf": "gtf", 
    //       "star_index": "star_index", 
    //       "extra_star_align_args": "extra_star_align_args", 
    //       "star_ignore_sjdbgtf": "star_ignore_sjdbgtf",
    //       "seq_platform": "seq_platform", 
    //       "seq_center": "seq_center",
    //       "with_umi": "with_umi", 
    //       "umi_dedup_stats": "umi_dedup_stats",
    //       "gtf_group_features": "gtf_group_features",
    //       "gtf_extra_attributes": "gtf_extra_attributes",
    //       "salmon_quant_libtype": "salmon_quant_libtype" 
    //     ],
    //     toState: [
    //       "star_alignment": "star_alignment", 
    //       "star_multiqc": "star_multiqc", 
    //       "genome_bam_indexed": "genome_bam_indexed", 
    //       "genome_bam_stats": "genome_bam_stats", 
    //       "genome_bam_flagstat": "genome_bam_flagstat", 
    //       "genome_bam_idxstats": "genome_bam_idxstats", 
    //       "transcriptome_bam_indexed": "transcriptome_bam_indexed", 
    //       "transcriptome_bam_stats": "transcriptome_bam_stats", 
    //       "transcriptome_bam_flagstat": "transcriptome_bam_flagstat", 
    //       "transcriptome_bam_idxstats": "transcriptome_bam_idxstats",
    //       "salmon_quant_results": "salmon_quant_output", 
    //       "salmon_tx2gene": "salmon_tx2gene_tsv", 
    //       "salmon_tximport": "salmon_tximport", 
    //       "salmon_summarizedexperiment": "salmon_summarizedexperiment" 
    //     ], 
    //     debug: true
    // )
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
