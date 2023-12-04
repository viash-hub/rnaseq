// Note: some helper functionality is added at the end of this file

workflow run_wf {
    take:
    input_ch

    main:
        output_ch = input_ch     

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
            fromState: [ 
                "input": "genome_bam", 
            ],
            toState: [ "genome_bam": "output" ],
            key: "genome_sorted_MarkDuplicates" 
        )
        | samtools_index.run (
            runIf: { id, state -> !state.skip_markduplicates && !state.with_umi },
            fromState: [
                "input": "genome_bam", 
                "bam_csi_index": "bam_csi_index" 
            ],
            toState: [
                "genome_bam_bai": "output_bai", 
                "genome_bam_csi": "output_csi" 
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
            "markduplicates_metrics": "markduplicates_metrics", 
            "stringtie_transcript_gtf": "stringtie_transcript_gtf",
            "stringtie_coverage_gtf": "stringtie_coverage_gtf",
            "stringtie_abundance": "stringtie_abundance",
            "stringtie_ballgown": "stringtie_ballgown", 
            // "featurecounts": "featurecounts",
            // "featurecounts_summary": "featurecounts_summary", 
            // "featurecounts_multiqc": "featurecounts_multiqc", 
            "bedgraph_forward": "bedgraph_forward",
            "bedgraph_reverse": "bedgraph_reverse",
            "bigwig_forward": "bigwig_forward",
            "bigwig_reverse": "bigwig_reverse"
        )

    emit:
        output_ch
}
