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
                "processed_genome_bam": "output_bam",
                "markduplicates_metrics": "metrics"
            ]
        )
        | samtools_sort.run (
            runIf: { id, state -> !state.skip_markduplicates && !state.with_umi },
            fromState: [ "input": "processed_genome_bam" ],
            toState: [ "processed_genome_bam": "output" ],
            key: "genome_sorted_MarkDuplicates" 
        )
        | samtools_index.run (
            runIf: { id, state -> !state.skip_markduplicates && !state.with_umi },
            fromState: [
                "input": "processed_genome_bam", 
                "csi": "bam_csi_index"
            ],
            toState: [ "genome_bam_index": "output" ],
            key: "genome_sorted_MarkDuplicates", 
        )
        | samtools_stats.run (
            runIf: { id, state -> !state.skip_markduplicates && !state.with_umi }, 
            fromState: [
                "input": "processed_genome_bam", 
                "bai": "genome_bam_index" 
            ],
            toState: [ "genome_bam_stats": "output" ],
            key: "MarkDuplicates_stats"
        )
        | samtools_flagstat.run (
            runIf: { id, state -> !state.skip_markduplicates && !state.with_umi }, 
            fromState: [
                "bam": "processed_genome_bam", 
                "bai": "genome_bam_index"
            ],
            toState: [ "genome_bam_flagstat": "output" ],
            key: "MarkDuplicates_flagstat"
        )
        | samtools_idxstats.run(
            runIf: { id, state -> !state.skip_markduplicates && !state.with_umi }, 
            fromState: [
            "bam": "processed_genome_bam", 
            "bai": "genome_bam_index"
            ],
            toState: [ "genome_bam_idxstats": "output" ],
            key: "MarkDuplicates_idxstats"
        ) 

        | stringtie.run (
            runIf: { id, state -> !state.skip_stringtie }, 
            fromState: [
                "strandedness": "strandedness", 
                "bam": "processed_genome_bam",
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
                "bam": "processed_genome_bam",
                "extra_bedtools_args": "extra_bedtools_args"
            ],
            toState: [
                "bedgraph_forward": "bedgraph_forward",
                "bedgraph_reverse": "bedgraph_reverse"
            ]
        )

        | bedclip.run (
            runIf: { id, state -> !state.skip_bigwig },
            fromState: [
                "input_bedgraph": "bedgraph_forward", 
                "sizes": "chrom_sizes" 
            ],
            toState: [ "bedgraph_forward": "output_bedgraph" ], 
            key: "bedclip_forward"
        )

        | bedgraphtobigwig.run (
            runIf: { id, state -> !state.skip_bigwig },
            fromState: [
                "bedgraph": "bedgraph_forward", 
                "sizes": "chrom_sizes" 
            ],
            toState: [ "bigwig_forward": "bigwig" ], 
            key: "bedgraphtobigwig_forward"
        )

        | bedclip.run (
            runIf: { id, state -> !state.skip_bigwig },
            fromState: [
                "input_bedgraph": "bedgraph_reverse", 
                "sizes": "chrom_sizes", 
            ],
            toState: [ "bedgraph_reverse": "output_bedgraph" ], 
            key: "bedclip_reverse"
        )

        | bedgraphtobigwig.run (
            runIf: { id, state -> !state.skip_bigwig },
            fromState: [
                "bedgraph": "bedgraph_reverse", 
                "sizes": "chrom_sizes" 
            ],
            toState: [ "bigwig_reverse": "bigwig" ], 
            key: "bedgraphtobigwig_reverse"
        )

        | map { id, state -> 
          def mod_state = state.findAll { key, value -> value instanceof java.nio.file.Path && value.exists() }
          [ id, mod_state ]
        }
        | niceView()

        | setState (
            "processed_genome_bam": "processed_genome_bam", 
            "genome_bam_index": "genome_bam_index",
            "genome_bam_stats": "genome_bam_stats",
            "genome_bam_flagstat": "genome_bam_flagstat", 
            "genome_bam_idxstats": "genome_bam_idxstats", 
            "markduplicates_metrics": "markduplicates_metrics", 
            "stringtie_transcript_gtf": "stringtie_transcript_gtf",
            "stringtie_coverage_gtf": "stringtie_coverage_gtf",
            "stringtie_abundance": "stringtie_abundance",
            "stringtie_ballgown": "stringtie_ballgown", 
            "bedgraph_forward": "bedgraph_forward",
            "bedgraph_reverse": "bedgraph_reverse",
            "bigwig_forward": "bigwig_forward",
            "bigwig_reverse": "bigwig_reverse"
        )

    emit:
        output_ch
}
