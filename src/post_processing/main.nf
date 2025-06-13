workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch     
      | picard_markduplicates.run (
        fromState: [
            "bam": "genome_bam",
            "fasta": "fasta",
            "fai": "fai",
        ], 
        toState: [
            "processed_genome_bam": "output_bam",
            "markduplicates_metrics": "metrics"
        ],
        directives: [ label: [ "midmem", "midcpu" ] ]
      )
      | samtools_sort.run (
        fromState: [ "input": "processed_genome_bam" ],
        toState: [ "processed_genome_bam": "output" ],
        key: "genome_sorted_MarkDuplicates",
        directives: [ label: [ "midmem", "midcpu" ] ]
      )
      | samtools_index.run (
        fromState: [
            "input": "processed_genome_bam", 
            "csi": "bam_csi_index"
        ],
        toState: [ "genome_bam_index": "output" ],
        key: "genome_sorted_MarkDuplicates",
        directives: [ label: [ "midmem", "midcpu" ] ]
      )
      | samtools_stats.run (
        fromState: [
            "input": "processed_genome_bam", 
            "bai": "genome_bam_index" 
        ],
        toState: [ "genome_bam_stats": "output" ],
        key: "MarkDuplicates_stats",
        directives: [ label: [ "midmem", "midcpu" ] ]
      )
      | samtools_flagstat.run (
        fromState: [
            "bam": "processed_genome_bam", 
            "bai": "genome_bam_index"
        ],
        toState: [ "genome_bam_flagstat": "output" ],
        key: "MarkDuplicates_flagstat",
        directives: [ label: [ "midmem", "midcpu" ] ]
      )
      | samtools_idxstats.run(
        fromState: [
          "bam": "processed_genome_bam", 
          "bai": "genome_bam_index"
        ],
        toState: [ "genome_bam_idxstats": "output" ],
        key: "MarkDuplicates_idxstats",
        directives: [ label: [ "midmem", "midcpu" ] ]
      ) 
      | setState (
        "processed_genome_bam": "processed_genome_bam", 
        "genome_bam_index": "genome_bam_index",
        "genome_bam_stats": "genome_bam_stats",
        "genome_bam_flagstat": "genome_bam_flagstat", 
        "genome_bam_idxstats": "genome_bam_idxstats", 
        "markduplicates_metrics": "markduplicates_metrics", 
      )

  emit:
      output_ch
}