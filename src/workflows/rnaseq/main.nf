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
    | map { id, state -> 
      (isBelowMaxContigSize(state.fai)) ? [id, state] : [id, state + [bam_csi_index: true]]
    }

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
          "umitools_bc_pattern": "umitools_bc_pattern",
          "umitools_bc_pattern2": "umitools_bc_pattern2",
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
          "fastqc_html_1": "fastqc_html_1",
          "fastqc_html_2": "fastqc_html_2",
          "fastqc_zip_1": "fastqc_zip_1",
          "fastqc_zip_2": "fastqc_zip_2",  
          "fastq_1": "qc_output1",
          "fastq_2": "qc_output2", 
          "trim_log_1": "trim_log_1", 
          "trim_log_2": "trim_log_2", 
          "trim_zip_1": "trim_zip_1",
          "trim_zip_2": "trim_zip_2",
          "trim_html_1": "trim_html_1",
          "trim_html_2": "trim_html_2",
          "passed_trimmed_reads": "passed_trimmed_reads",
          "num_trimmed_reads": "num_trimmed_reads",
          "sortmerna_log": "sortmerna_log",
          "salmon_json_info": "salmon_json_info"
        ]
    )

    // Infer strandedness from Salmon pseudo-alignment results
    | map { id, state -> 
    (state.strandedness == 'auto') ? 
      [ id, state + [strandedness: getSalmonInferredStrandedness(state.salmon_json_info)] ] : 
      [id, state] 
    }

    // Filter FastQ files based on minimum trimmed read count after adapter trimming
    | map { id, state -> 
      def input = state.fastq_2 ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
      def num_reads = state.min_trimmed_reads + 1
      num_reads = 
        (!state.skip_trimming && input.size() == 2) ?
          getTrimGaloreReadsAfterFiltering(state.trim_log_2) : 
          getTrimGaloreReadsAfterFiltering(state.trim_log_1)
      def passed_trimmed_reads = 
        (state.skip_trimming || (num_reads >= state.min_trimmed_reads)) ? 
          true : 
          false 
      [ id, state + [num_trimmed_reads: num_reads, passed_trimmed_reads: passed_trimmed_reads] ] 
    }
    // | filter { id, state -> state.skip_trimming || state.passed_trimmed_reads }
    // TODO: Get list of samples that failed trimming threshold for MultiQC report

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
          "salmon_quant_results": "salmon_quant_results"
        ]
    )

    // Filter channels to get samples that passed STAR minimum mapping percentage
    | map { id, state -> 
      def percent_mapped = getStarPercentMapped(state.star_multiqc) 
      def passed_mapping = (percent_mapped >= state.min_mapped_reads) ? true : false
      [ id, state + [percent_mapped: percent_mapped, passed_mapping: passed_mapping] ]
    }
    // | filter { id, state -> state.passed_mapping) }
    // TODO: Get list of samples that failed mapping for MultiQC report
    
    // Post-processing
    // | post_processing.run (
    //     fromState: [
    //       "id": "id", 
    //       "paired": "paired", 
    //       "strandedness": "strandedness", 
    //       "fasta": "fasta",
    //       "fai": "fai", 
    //       "gtf": "gtf", 
    //       "genome_bam": "genome_bam_sorted", 
    //       "chrom_sizes": "chrom_sizes", 
    //       "star_multiqc": "star_multiqc",
    //       "extra_picard_args": "extra_picard_args", 
    //       "extra_stringtie_args": "extra_stringtie_args", 
    //       "stringtie_ignore_gtf": "stringtie_ignore_gtf", 
    //       "extra_bedtools_args": "extra_bedtools_args", 
    //       "extra_featurecounts_args": "extra_featurecounts_args", 
    //       "bam_csi_index": "bam_csi_index", 
    //       "min_mapped_reads": "min_mapped_reads", 
    //       "with_umi": "with_umi",
    //       "biotype": "biotype", 
    //       "biotypes_header": "biotypes_header",
    //       "skip_qc": "skip_qc",
    //       "skip_markdupkicates": "skip_markdupkicates", 
    //       "skip_stringtie": "skip_stringtie", 
    //       "skip_biotype_qc": "skip_biotype_qc", 
    //       "skip_bigwig": "skip_bigwig", 
    //       "featurecounts_group_type": "featurecounts_group_type", 
    //       "featurecounts_feature_type": "featurecounts_feature_type", 
    //       "gencode": "gencode"
    //     ], 
    //     toState: [
    //       "genome_bam_sorted": "processed_genome_bam", 
    //       "genome_bam_index": "genome_bam_index",
    //       "genome_bam_stats": "genome_bam_stats",
    //       "genome_bam_flagstat": "genome_bam_flagstat", 
    //       "genome_bam_idxstats": "genome_bam_idxstats", 
    //       "markduplicates_metrics": "markduplicates_metrics",
    //       "stringtie_transcript_gtf": "stringtie_transcript_gtf",
    //       "stringtie_coverage_gtf": "stringtie_coverage_gtf",
    //       "stringtie_abundance": "stringtie_abundance",
    //       "stringtie_ballgown": "stringtie_ballgown", 
    //       "featurecounts": "featurecounts",
    //       "featurecounts_summary": "featurecounts_summary", 
    //       "featurecounts_multiqc": "featurecounts_multiqc", 
    //       "bedgraph_forward": "bedgraph_forward",
    //       "bedgraph_reverse": "bedgraph_reverse",
    //       "bigwig_forward": "bigwig_forward",
    //       "bigwig_reverse": "bigwig_reverse"
    //     ], 
    // )

    | salmon_quant_merge_counts.run (
        fromState: [ 
          "salmon_quant_results": "salmon_quant_results", 
          "gtf": "gtf", 
          "gtf_extra_attributes": "gtf_extra_attributes", 
          "gtf_group_features": "gtf_group_features"],
        toState: [
          "tpm_gene": "tpm_gene",
          "counts_gene": "counts_gene",
          "counts_gene_length_scaled": "counts_gene_length_scaled",
          "counts_gene_scaled": "counts_gene_scaled", 
          "tpm_transcript": "tpm_transcript", 
          "counts_transcript": "counts_transcript", 
          "salmon_merged_summarizedexperiment": "salmon_merged_summarizedexperiment"
        ]
    )

    // Final QC
    // | final_qc.run (
    //     fromState: [
    //       "id": "id", 
    //       "paired": "paired", 
    //       "strandedness": "strandedness", 
    //       "gtf": "gtf", 
    //       "genome_bam": "genome_bam_sorted", 
    //       "genome_bam_index": "genome_bam_index",
    //       "salmon_quant_results": "salmon_quant_output", 
    //       "gene_bed": "gene_bed",
    //       "extra_preseq_args": "extra_preseq_args", 
    //       "extra_deseq2_args": "extra_deseq2_args",
    //       "extra_deseq2_args2": "extra_deseq2_args2",
    //     ], 
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
def isBelowMaxContigSize(fai_file) {
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
      return false
    }
  }
  return true
}

import groovy.json.JsonSlurper
//
// Function that parses Salmon quant 'meta_info.json' output file to get inferred strandedness
//
def getSalmonInferredStrandedness(json_file) {
  def lib_type = new JsonSlurper().parseText(json_file.text).get('library_types')[0]
  def strandedness = 'reverse'
  if (lib_type) {
    if (lib_type in ['U', 'IU']) {
      strandedness = 'unstranded'
    } 
    else if (lib_type in ['SF', 'ISF']) {
      strandedness = 'forward'
    } 
    else if (lib_type in ['SR', 'ISR']) {
      strandedness = 'reverse'
    }
  }
  return strandedness
}

//
// Function that parses TrimGalore log output file to get total number of reads after trimming
//
def getTrimGaloreReadsAfterFiltering(log_file) {
  def total_reads = 0
  def filtered_reads = 0
  log_file.eachLine { line ->
    def total_reads_matcher = line =~ /([\d\.]+)\ssequences processed in total/
    def filtered_reads_matcher = line =~ /shorter than the length cutoff[^:]+:\s([\d\.]+)/
    if (total_reads_matcher) total_reads = total_reads_matcher[0][1].toFloat()
    if (filtered_reads_matcher) filtered_reads = filtered_reads_matcher[0][1].toFloat()
  }
  return total_reads - filtered_reads
}

//
// Function that parses and returns the alignment rate from the STAR log output
//
def getStarPercentMapped(align_log) {
  def percent_aligned = 0
  def pattern = /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/
  align_log.eachLine { line ->
    def matcher = line =~ pattern
    if (matcher) {
        percent_aligned = matcher[0][1].toFloat()
    }
  }
  return percent_aligned
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