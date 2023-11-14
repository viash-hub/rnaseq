// Note: some helper functionality is added at the end of this file

workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch   
    
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
        key: "genome_sorted_MarkDuplicates" 
    )
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

    // Feature biotype QC using featureCounts
    // PREPARE_GENOME
    //         .out
    //         .gtf
    //         .map { WorkflowRnaseq.biotypeInGtf(it, biotype, log) }
    //         .set { biotype_in_gtf }

    //     // Prevent any samples from running if GTF file doesn't have a valid biotype
    //     ch_genome_bam
    //         .combine(PREPARE_GENOME.out.gtf)
    //         .combine(biotype_in_gtf)
    //         .filter { it[-1] }
    //         .map { it[0..<it.size()-1] }
    //         .set { ch_featurecounts }

    | subread_featurecounts.run (
        runIf: { id, state -> !state.skip_qc && !state.skip_biotype_qc && state.biotype },
        fromState: [
            "paired": "paired", 
            "strandedness": "strandedness", 
            "gtf": "gtf", 
            "bam": "genome_bam", 
            "gencode": "gencode",
            "extra_featurecounts_args": "extra_featurecounts_args",
            "featurecounts_group_type": "featurecounts_group_type",
            "featruecount_feature_type": "featruecount_feature_type",
        ],
        toState: [
            "featurecounts": "counts",
            "featurecounts_summary": "summary"
        ]
    )

    | multiqc_custom_biotype.run (
        runIf: { id, state -> !state.skip_qc && !state.skip_biotype_qc && state.biotype },
        fromState: [
            "id": "id",
            "biocounts": "featurecounts", 
            "biotypes_header": "biotypes_header"
        ],
        toState: [ "featurecounts_multiqc": "featurecounts_multiqc" ]
    )

    // Genome-wide coverage with BEDTools

    | bedtools_genomecov.run (
        runIf: { id, state -> !state.skip_bigwig },
        fromState: [
            "strandedness": "strandedness", 
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
        toState: [ "bigwig_forward": "bigwig" ], 
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
    )

    | view { "Output: $it" }

  emit:
    output_ch
}

//
// Function to check whether biotype field exists in GTF file
//
def biotypeInGtf(gtf_file, biotype, log) {
    def hits = 0
    gtf_file.eachLine { line ->
        def attributes = line.split('\t')[-1].split()
        if (attributes.contains(biotype)) {
            hits += 1
        }
    }
    if (hits) {
        return true
    } else {
        log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Biotype attribute '${biotype}' not found in the last column of the GTF file!\n\n" +
            "  Biotype QC will be skipped to circumvent the issue below:\n" +
            "  https://github.com/nf-core/rnaseq/issues/460\n\n" +
            "  Amend '--featurecounts_group_type' to change this behaviour.\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        return false
    }
}