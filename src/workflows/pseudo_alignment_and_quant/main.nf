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

    // Count reads from BAM alignments using Salmon
    | salmon_quant.run ( 
        runIf: { id, state -> state.pseudo_aligner == 'salmon'}
        fromState: [
          "input": "input", 
          "transcript_fasta": "transcript_fasta", 
          "gtf": "gtf", 
          "salmon_index": "pseudo_index",
          "versions": "versions" 
        ],
        args: ["alignment_mode": false],
        toState: [
          "quant_results": "output", 
          "pseudo_multiqc": "output"
          "versions": "updated_versions"
        ]
    )

    | kallisto_quant.run ( 
        runIf: { id, state -> state.pseudo_aligner == 'kallisto'}
        fromState: [
          "input": "input", 
          "gtf": "gtf", 
          "index": "pseudo_index",
          "fragment_length": "kallisto_quant_fragment_length", 
          "fragment_length_sd": "kallisto_quant_fragment_length_sd", 
          "versions": "versions" 
        ],
        toState: [
          "quant_results": "output", 
          "pseudo_multiqc": "log"
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
        "quant_results": "quant_results", 
        "updated_versions": "versions" ]
    )

  emit:
    output_ch
}
