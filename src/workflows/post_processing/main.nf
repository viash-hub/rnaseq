// Note: some helper functionality is added at the end of this file

workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch

    // TODO: Move to final_qc
    // | deseq2_qc.run (
    //     runIf: { id, state -> !state.skip_qc && !skip_deseq2_qc },
    //     fromState: [
    //         "counts": "counts_gene_length_scaled",
    //         "pca_header_multiqc": "pca_header_multiqc", 
    //         "clustering_header_multiqc": "clustering_header_multiqc",
    //         "deseq2_vst": "deseq2_vst"
    //     ], 
    //     toState: [
    //         "deseq2_output": "deseq2_output", 
    //         "deseq2_pca_multiqc": "pca_multiqc", 
    //         "deseq2_dists_multiqc": "dists_multiqc"
    //     ]
    // )
    
    // Filter channels to get samples that passed STAR minimum mapping percentage
    // | map_star_pass_fail.run (
    //     fromState: [
    //         "star_multiqc": "star_multiqc", 
    //         "min_mapped_reads": "min_mapped_reads", 
    //         "genome_bam": "genome_bam"
    //     ], 
    //     toState: [
    //         "percent_mapped": "percent_mapped",
    //         "star_passed": "star_pass"
    //     ]
    // )
    // | filter { id, state -> state.star_passed }
    // TODO: star_failed_multiqc: multiqcTsvFromList

    // TODO: Move to final_qc
    // | preseq_lcextrap.run (
    //     runIf: { id, state -> !state.skip_qc && !skip_preseq },
    //     fromState: [
    //         "paired": "paired",
    //         "bam": "genome_bam",
    //         "extra_preseq_arg": "extra_preseq_arg"
    //     ],
    //     toState: [ "preseq_output": "output" ],
    // )
    
    | picard_markduplicates.run (
        runIf: { id, state -> !state.skip_markduplicates && !state.with_umi }, 
        fromState: [
            "bam": "genome_bam",
            "fasta": "fasta",
            "fai": "fai",
            "extra_picard_args": "extra_picard_args"
        ], 
        toState: [
            "genome_bam": "output_bam",
            "markduplicates_metrics": "metrics"
        ]
    )
    | samtools_sort.run (
        runIf: { id, state -> !state.skip_markduplicates && !state.with_umi },
        fromState: ["input": "genome_bam"],
        toState: ["genome_bam": "output"],
        key: "genome_sorted_MarkDuplicates"    )
    | samtools_index.run (
      runIf: { id, state -> !state.skip_markduplicates && !state.with_umi },
      fromState: [
        "bam": "genome_bam", 
        "bam_csi_index": "bam_csi_index" 
      ],
      toState: [
        "genome_bam_bai": "bai", 
        "genome_bam_csi": "csi" 
      ],
      key: "genome_sorted_MarkDuplicates", 
    )
    | map { id, state -> 
      def genome_bam_index = state.genome_bam_bai ? state.genome_bam_bai : state.genome_bam_csi
      [ id, state + [ genome_bam_index: genome_bam_index ] ]
    }
    | samtools_stats.run (
        runIf: { id, state -> !state.skip_markduplicates && !state.with_umi }, 
        fromState: [
          "bam": "genome_bam", 
          "bai": "genome_bam_index" 
        ],
        toState: ["genome_bam_stats": "output"],
        key: "MarkDuplicates_stats"
    )
    | samtools_flagstat.run (
        runIf: { id, state -> !state.skip_markduplicates && !state.with_umi }, 
        fromState: [
          "bam": "genome_bam", 
          "bai": "genome_bam_index" 
        ],
        toState: ["genome_bam_flagstat": "output"],
        key: "MarkDuplicates_flagstat"
    )
    | samtools_idxstats.run(
        runIf: { id, state -> !state.skip_markduplicates && !state.with_umi }, 
        fromState: [
          "bam": "genome_bam", 
          "bai": "genome_bam_index" 
        ],
        toState: ["genome_bam_idxstats": "output"],
        key: "MarkDuplicates_idxstats"
    ) 

    | stringtie.run (
        runIf: { id, state -> !state.skip_stringtie }, 
        fromState: [
            "strandedness": "strandedness", 
            "bam": "genome_bam",
            "annotation_gtf": "gtf",
            "extra_stringtie_args": "extra_stringtie_args"
        ], 
        toState: [
            "stringtie_transcript_gtf": "transcript_gtf",
            "stringtie_coverage_gtf": "coverage_gtf",
            "stringtie_abundance": "abundance",
            "stringtie_ballgown": "ballgown"
        ]
    )

    // // Feature biotype QC using featureCounts
    // | check_biotype_in_gtf.run (
    //   runIf: { id, state -> !state.skip_qc && !state.skip_biotype_qc && state.biotype },
    //     fromState: [
    //         "gtf": "gtf",
    //         "biotype": "biotype"
    //     ],
    //     toState: [ "biotype_exists": "biotype_exists" ]
    // )

    // | subread_featurecounts.run (
    //     runIf: { id, state -> !state.skip_qc && !state.skip_biotype_qc && state.biotype },
    //     fromState: [
    //         "paired": "paired", 
    //         "strandedness": "strandedness", 
    //         "gtf": "gtf", 
    //         "bam": "genome_bam"
    //     ],
    //     toState: [
    //         "featurecounts": "counts",
    //         "featurecounts_summary": "summary"
    //     ]
    // )

    // | multiqc_custom_biotype.run (
    //     runIf: { id, state -> !state.skip_qc && !state.skip_biotype_qc && state.biotype },
    //     fromState: [
    //         "id": "id",
    //         "biocounts": "featurecounts", 
    //         "biotypes_header": "biotypes_header"
    //     ],
    //     toState: [ "featurecounts_multiqc": "featurecounts_multiqc" ]
    // )

    // Genome-wide coverage with BEDTools

    | bedtools_genomecov.run (
        runIf: { id, state -> !state.skip_bigwig },
        fromState: [
            "bam": "genome_bam",
            "extra_bedtools_args": "extra_bedtools_args"
        ],
        toState: [
            "bedgraph_forward": "bedgraph_forward",
            "bedgraph_reverse": "bedgraph_reverse"
        ]
    )

    | ucsc_bedclip.run (
        runIf: { id, state -> !state.skip_bigwig },
        fromState: [
            "input_bedgraph": "bedgraph_forward", 
            "sizes": "chrom_sizes"
        ],
        toState: [ "bedgraph_forward": "output_bedgraph" ], 
        key: "bedclip_forward"
    )

    | ucsc_bedgraphtobigwig.run (
        runIf: { id, state -> !state.skip_bigwig },
        fromState: [
            "bedgraph": "bedgraph_forward", 
            "sizes": "chrom_sizes"
        ],
        toState: [ " bigwig_forward": "bigwig" ], 
        key: "bedgraphtobigwig_forward"
    )

    | ucsc_bedclip.run (
        runIf: { id, state -> !state.skip_bigwig },
        fromState: [
            "input_bedgraph": "bedgraph_reverse", 
            "sizes": "chrom_sizes"
        ],
        toState: [ "bedgraph_reverse": "output_bedgraph" ], 
        key: "bedclip_reverse"
    )

    | ucsc_bedgraphtobigwig.run (
        runIf: { id, state -> !state.skip_bigwig },
        fromState: [
            "bedgraph": "bedgraph_reverse", 
            "sizes": "chrom_sizes"
        ],
        toState: [ "bigwig_reverse": "bigwig" ], 
        key: "bedgraphtobigwig_reverse"
    )

    | setState (
        "processed_genome_bam": "genome_bam", 
        "genome_bam_indexed": "genome_bam_indexed",
        "genome_bam_stats": "genome_bam_stats",
        "genome_bam_flagstat": "genome_bam_flagstat", 
        "genome_bam_idxstats": "genome_bam_idxstats", 
        "stringtie_transcript_gtf": "stringtie_transcript_gtf",
        "transcript_gtf": "transcript_gtf",
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
    )

    | view { "Output: $it" }

  emit:
    output_ch
}

