workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch

    | map { id, state ->
      def biotype = state.gencode ? "gene_type" : state.featurecounts_group_type 
      def bc_pattern = state.umitools_bc_pattern2 ? [ state.umitools_bc_pattern, state.umitools_bc_pattern2 ] : [ state.umitools_bc_pattern ] 
      [ id, state + [ biotype: biotype, bc_pattern: bc_pattern ] ]
    } 

    // prepare all the necessary files for reference genome
    | prepare_genome.run ( 
        fromState: [
          "fasta": "fasta", 
          "gtf": "gtf", 
          "gff": "gff",
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
        ]
    )

    // Check if contigs in genome fasta file > 512 Mbp
    // | map { checkMaxContigSize(it[1].fai) }
    // | niceView()

    // Concatenate FastQ files from same sample if required
    | cat_fastq.run (
        fromState: [
          "read_1": "fastq_1", 
          "read_2": "fastq_2"
        ], 
        toState: [ 
          "fastq_1": "fastq_1",
          "fastq_2": "fastq_2"
        ]
    )
    
    // Pre-process fastq files
    | pre_processing.run ( 
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
          "ribo_database_manifest": "ribo_database_manifest", 
          "salmon_index": "salmon_index"
        ], 
        toState: [ 
          "fastqc_report": "fastqc_report", 
          "fastq_1": "qc_output1",
          "fastq_2": "qc_output2", 
          "trim_log": "trim_log", 
          "trim_zip": "trim_zip",
          "trim_html": "trim_html",
          "sortmerna_log": "sortmerna_log"
        ]
    )

    // Genome alignment and quantification
    | genome_alignment_and_quant.run (
        fromState: [
          "id": "id", 
          "fastq_1": "fastq_1",
          "fastq_2": "fastq_2", 
          "strandedness": "strandedness", 
          "gtf": "gtf",
          "transcript_fasta": "transcript_fasta",
          "bam_csi_index": "bam_csi_index", 
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
          "genome_bam_sorted": "genome_bam_sorted",
          "genome_bam_index": "genome_bam_index", 
          "genome_bam_stats": "genome_bam_stats", 
          "genome_bam_flagstat": "genome_bam_flagstat", 
          "genome_bam_idxstats": "genome_bam_idxstats", 
          "transcriptome_bam_sorted": "transcriptome_bam_sorted", 
          "transcriptome_bam_index": "transcriptome_bam_index", 
          "transcriptome_bam_stats": "transcriptome_bam_stats", 
          "transcriptome_bam_flagstat": "transcriptome_bam_flagstat", 
          "transcriptome_bam_idxstats": "transcriptome_bam_idxstats",
          "salmon_quant_results": "salmon_quant_output"
        ]
    )
    
    // Post-processing
    | post_processing.run (
        fromState: [
          "id": "id", 
          "paired": "paired", 
          "strandedness": "strandedness", 
          "fasta": "fasta",
          "fai": "fai", 
          "gtf": "gtf", 
          "genome_bam": "genome_bam_sorted", 
          "chrom_sizes": "chrom_sizes", 
          "star_multiqc": "star_multiqc",
          "extra_picard_args": "extra_picard_args", 
          "extra_stringtie_args": "extra_stringtie_args", 
          "stringtie_ignore_gtf": "stringtie_ignore_gtf", 
          "extra_bedtools_args": "extra_bedtools_args", 
          "extra_featurecounts_args": "extra_featurecounts_args", 
          "bam_csi_index": "bam_csi_index", 
          "min_mapped_reads": "min_mapped_reads", 
          "with_umi": "with_umi",
          "biotype": "biotype", 
          "biotypes_header": "biotypes_header",
          "skip_qc": "skip_qc",
          "skip_markdupkicates": "skip_markdupkicates", 
          "skip_stringtie": "skip_stringtie", 
          "skip_biotype_qc": "skip_biotype_qc", 
          "skip_bigwig": "skip_bigwig", 
          "featurecounts_group_type": "featurecounts_group_type", 
          "featurecounts_feature_type": "featurecounts_feature_type", 
          "gencode": "gencode"
        ], 
        toState: [
          "genome_bam_sorted": "processed_genome_bam", 
          "genome_bam_index": "genome_bam_index",
          "genome_bam_stats": "genome_bam_stats",
          "genome_bam_flagstat": "genome_bam_flagstat", 
          "genome_bam_idxstats": "genome_bam_idxstats", 
          "stringtie_transcript_gtf": "stringtie_transcript_gtf",
          "stringtie_coverage_gtf": "stringtie_coverage_gtf",
          "stringtie_abundance": "stringtie_abundance",
          "stringtie_ballgown": "stringtie_ballgown", 
          "featurecounts": "featurecounts",
          "featurecounts_summary": "featurecounts_summary", 
          "featurecounts_multiqc": "featurecounts_multiqc", 
          "bedgraph_forward": "bedgraph_forward",
          "bedgraph_reverse": "bedgraph_reverse",
          "bigwig_forward": "bigwig_forward",
          "bigwig_reverse": "bigwig_reverse"
        ], 
    )

    // Final QC
    // | final_qc.run (
    //     fromState: [], 
    //     toState: [] 
    // )

    | niceView()

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

import nextflow.Nextflow

//
// Function to generate an error if contigs in genome fasta file > 512 Mbp
//
def checkMaxContigSize(fai_file) {
  def max_size = 512000000
  fai_file.eachLine { line ->
    def lspl  = line.split('\t')
    def chrom = lspl[0]
    def size  = lspl[1]
    if (size.toInteger() > max_size) {
      def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
          "  Contig longer than ${max_size}bp found in reference genome!\n\n" +
          "  ${chrom}: ${size}\n\n" +
          "  Provide the '--bam_csi_index' parameter to use a CSI instead of BAI index.\n\n" +
          "  Please see:\n" +
          "  https://github.com/nf-core/rnaseq/issues/744\n" +
          "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
      Nextflow.error(error_string)
    }
  }
}