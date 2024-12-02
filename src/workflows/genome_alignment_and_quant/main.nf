// TODO: Improve logic of this wf, e.g. by splitting up in 2 subwfs
// TODO: See if this can be aligned with the pseudo-alignment branch of the logic
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

      | star_align_reads.run(
        runIf: { id, state -> state.aligner == 'star_salmon' },
        fromState: [
          "input": "fastq_1",
          "input_r2": "fastq_2",
          "genome_dir": "star_index",
          "sjdb_gtf_file": "gtf",
          "out_sam_attr_rg_line": "star_sam_attr_rg_line",
          "sjdb_gtf_feature_exon": "star_sjdb_gtf_feature_exon"
        ],
        toState: [
          "genome_bam": "aligned_reads",
          "transcriptome_bam": "reads_aligned_to_transcriptome",
          "star_multiqc": "log"
        ],
        args: [ 
          quant_mode: "TranscriptomeSAM", 
          twopass_mode: "Basic", 
          out_sam_type: "BAM;Unsorted", 
          run_rng_seed: 0, 
          out_filter_multimap_nmax: 20, 
          align_sjdb_overhang_min: 1, 
          out_sam_attributes: "NH;HI;AS;NM;MD", 
          quant_transcriptome_sam_output: "BanSingleEnd" 
        ],
        directives:
          - label: [ highmem, midcpu ]
      )

      // GENOME BAM
      | samtools_sort.run(
        runIf: { id, state -> state.aligner == 'star_salmon' },
        fromState: ["input": "genome_bam"],
        toState: ["genome_bam_sorted": "output"],
        key: "genome_sorted"
        directives:
          - label: [ midmem, midcpu ]
      )
      | samtools_index.run(
        runIf: { id, state -> state.aligner == 'star_salmon' },
        fromState: [ 
          "input": "genome_bam_sorted", 
          "csi": "bam_csi_index"
        ],
        toState: [ "genome_bam_index": "output" ],
        key: "genome_sorted"
        directives:
          - label: [ midmem, midcpu ]
      )
      | samtools_stats.run(
        runIf: { id, state -> state.aligner == 'star_salmon' },
        fromState: [
          "input": "genome_bam_sorted", 
          "bai": "genome_bam_index", 
          "fasta": "fasta"
        ],
        toState: [ "genome_bam_stats": "output" ],
        key: "genome_stats"
        directives:
          - label: [ midmem, midcpu ]
      )
      | samtools_flagstat.run(
        runIf: { id, state -> state.aligner == 'star_salmon' },
        fromState: [
          "bam": "genome_bam_sorted", 
          "bai": "genome_bam_index", 
          "fasta": "fasta"
        ],
        toState: [ "genome_bam_flagstat": "output" ],
        key: "genome_flagstat"
        directives:
          - label: [ midmem, midcpu ]
      )
      | samtools_idxstats.run(
        runIf: { id, state -> state.aligner == 'star_salmon' },
        fromState: [
          "bam": "genome_bam_sorted", 
          "bai": "genome_bam_index", 
          "fasta": "fasta"
        ],
        toState: [ "genome_bam_idxstats": "output" ],
        key: "genome_idxstats"
        directives:
          - label: [ midmem, midcpu ]
      )

      //
      // Remove duplicate reads from BAM file based on UMIs
      // 
      
      // Deduplicate genome BAM file
      | umi_tools_dedup.run( 
        runIf: { id, state -> state.with_umi && state.aligner == 'star_salmon' },
        fromState: { id, state -> 
          def output_stats = state.umi_dedup_stats ? state.id : 
          [ paired: state.paired, 
          input: state.genome_bam, 
          bai: state.genome_bam_index,
          output_stats: output_stats] 
        },
        toState: [ "genome_bam_sorted": "output" ],
        key: "genome_deduped"
        directives:
          - label: [ midmem, midcpu ]
      )
      | samtools_index.run(
        runIf: { id, state -> state.with_umi && state.aligner == 'star_salmon' },
        fromState: [
          "input": "genome_bam_sorted", 
          "csi": "bam_csi_index"
        ],
        toState: [ "genome_bam_index": "output" ],
        key: "genome_deduped"
        directives:
          - label: [ midmem, midcpu ]
      )
      | samtools_stats.run(
        runIf: { id, state -> state.with_umi && state.aligner == 'star_salmon' },
        fromState: [
          "input": "genome_bam_sorted", 
          "bai": "genome_bam_index", 
          "fasta": "fasta"
        ],
        toState: [ "genome_bam_stats": "output" ],
        key: "genome_deduped_stats"
        directives:
          - label: [ midmem, midcpu ]
      )
      | samtools_flagstat.run(
        runIf: { id, state -> state.with_umi && state.aligner == 'star_salmon' },
        fromState: [
          "bam": "genome_bam_sorted", 
          "bai": "genome_bam_index", 
          "fasta": "fasta"
        ],
        toState: [ "genome_bam_flagstat": "output" ],
        key: "genome_deduped_flagstat"
        directives:
          - label: [ midmem, midcpu ]
      )
      | samtools_idxstats.run(
        runIf: { id, state -> state.with_umi && state.aligner == 'star_salmon' },
        fromState: [
          "bam": "genome_bam_sorted", 
          "bai": "genome_bam_index", 
          "fasta": "fasta", 
        ],
        toState: [ "genome_bam_idxstats": "output" ],
        key: "genome_deduped_idxstats"
        directives:
          - label: [ midmem, midcpu ]
      )

      // Deduplicate transcriptome BAM file

      | samtools_sort.run(
        runIf: { id, state -> state.with_umi && state.aligner == 'star_salmon' },
        fromState: [ "input": "transcriptome_bam" ],
        toState: [ "transcriptome_bam": "output" ],
        key: "transcriptome_sorted"
        directives:
          - label: [ midmem, midcpu ]
      )
      | samtools_index.run(
        runIf: { id, state -> state.with_umi && state.aligner == 'star_salmon' },
        fromState: [
          "input": "transcriptome_bam", 
          "csi": "bam_csi_index"
        ],
        toState: [ "transcriptome_bam_index": "output" ],
        key: "transcriptome_sorted"
        directives:
          - label: [ midmem, midcpu ]
      )
      | samtools_stats.run(
        runIf: { id, state -> state.with_umi && state.aligner == 'star_salmon' },
        fromState: [
          "input": "transcriptome_bam", 
          "bai": "transcriptome_bam_index", 
        ],
        toState: [ "transcriptome_bam_stats": "output" ],
        key: "transcriptome_stats"
        directives:
          - label: [ midmem, midcpu ]
      )
      | samtools_flagstat.run(
        runIf: { id, state -> state.with_umi && state.aligner == 'star_salmon' },
        fromState: [
          "bam": "transcriptome_bam", 
          "bai": "transcriptome_bam_index"
        ],
        toState: [ "transcriptome_bam_flagstat": "output" ],
        key: "transcriptome_flagstat"
        directives:
          - label: [ midmem, midcpu ]
      )
      | samtools_idxstats.run(
        runIf: { id, state -> state.with_umi && state.aligner == 'star_salmon' },
        fromState: [
          "bam": "transcriptome_bam", 
          "bai": "transcriptome_bam_index" 
        ],
        toState: [ "transcriptome_bam_idxstats": "output" ],
        key: "transcriptome_idxstats"
        directives:
          - label: [ midmem, midcpu ]
      )    
      
      | umi_tools_dedup.run( 
        runIf: { id, state -> state.with_umi && state.aligner == 'star_salmon' },
        fromState: { id, state -> 
          def output_stats = state.umi_dedup_stats ? state.id : 
          [ paired: state.paired, 
          input: state.transcriptome_bam, 
          bai: state.transcriptome_bam_index,
          output_stats: output_stats] 
        },
        toState: [ "transcriptome_bam_deduped": "output" ],
        key: "transcriptome_deduped"
        directives:
          - label: [ midmem, midcpu ]
      )
      | samtools_sort.run(
        runIf: { id, state -> state.with_umi && state.aligner == 'star_salmon' }, 
        fromState: [ "input": "transcriptome_bam_deduped" ],
        toState: [ "transcriptome_bam": "output" ],
        key: "transcriptome_deduped_sorted"
        directives:
          - label: [ midmem, midcpu ]
      )
      | samtools_index.run(
        runIf: { id, state -> state.with_umi && state.aligner == 'star_salmon' },
        fromState: [
          "input": "transcriptome_bam", 
          "csi": "bam_csi_index"
        ],
        toState: [ "transcriptome_bam_index": "output" ],
        key: "transcriptome_deduped_sorted"
        directives:
          - label: [ midmem, midcpu ]
      )
      | samtools_stats.run(
        runIf: { id, state -> state.with_umi && state.aligner == 'star_salmon' }, 
        fromState: [
          "input": "transcriptome_bam", 
          "bai": "transcriptome_bam_index"
        ],
        toState: [ "transcriptome_bam_stats": "output" ],
        key: "transcriptome_deduped_stats"
        directives:
          - label: [ midmem, midcpu ]
      )
      | samtools_flagstat.run(
        runIf: { id, state -> state.with_umi && state.aligner == 'star_salmon' }, 
        fromState: [
          "bam": "transcriptome_bam", 
          "bai": "transcriptome_bam_index"
        ],
        toState: [ "transcriptome_bam_flagstat": "output" ],
        key: "transcriptome_deduped_flagstat"
        directives:
          - label: [ midmem, midcpu ]
      )
      | samtools_idxstats.run(
        runIf: { id, state -> state.with_umi && state.aligner == 'star_salmon' }, 
        fromState: [
          "bam": "transcriptome_bam", 
          "bai": "transcriptome_bam_index" 
        ],
        toState: [ "transcriptome_bam_idxstats": "output" ],
        key: "transcriptome_deduped_idxstats"
        directives:
          - label: [ midmem, midcpu ]
      ) 

      // Fix paired-end reads in name sorted BAM file
      | umi_tools_prepareforrsem.run(
        runIf: { id, state -> state.with_umi && state.paired && state.aligner == 'star_salmon' },
        fromState: [ "input": "transcriptome_bam" ],
        toState: [ "transcriptome_bam": "output" ]
        directives:
          - label: [ midmem, midcpu ]
      )

      // Infer lib-type for salmon quant
      | map { id, state -> 
        def lib_type = (state.paired) ? 
          (
            (state.strandedness == "forward") ? 
              "ISF" : 
              (
                (state.strandedness == "reverse") ? "ISR" : "IU"
              )
          ) 
          : (
            (state.strandedness == "forward") ? 
              "SF" : 
              (
                (state.strandedness == "reverse") ? "SR" : "U"
              )
          ) 
        [ id, state + [lib_type: lib_type] ]
      }

      // Count reads from BAM alignments using Salmon
      | salmon_quant.run( 
        runIf: { id, state -> state.aligner == 'star_salmon' },
        fromState: [
          "lib_type": "lib_type",
          "alignments": "transcriptome_bam", 
          "targets": "transcript_fasta", 
          "gene_map": "gtf"
        ],
        toState: [ 
          "quant_out_dir": "output",
          "quant_results_file": "quant_results" 
        ]
        directives:
          - label: [ midmem, midcpu ]
      )

      | map { id, state -> 
        def mod_state = (state.aligner == 'star_salmon') ? state + [salmon_multiqc: state.quant_out_dir] : state
        [ id, mod_state ]
      }
    
      | rsem_calculate_expression.run(
        runIf: { id, state -> state.aligner == 'star_rsem' },
        fromState: [
          "id": "id",
          "strandedness": "strandedness",
          "paired": "paired",
          "input": "input",
          "index": "rsem_index", 
          "counts_gene": "rsem_counts_gene",
          "counts_transcripts": "rsem_counts_transcripts",
          "stat": "rsem_multiqc",
          "logs": "star_multiqc",
          "bam_star": "bam_star_rsem",
          "bam_genome": "bam_genome_rsem",
          "bam_transcript": "bam_transcript_rsem"          
        ],
        args: [
          star: true,
          star_output_genome_bam: true,
          star_gzipped_read_file: true,
          estimate_rspd: true,
          seed: 1
        ],
        toState: [
          "rsem_counts_gene": "counts_gene",
          "rsem_counts_transcripts": "counts_transcripts",
          "rsem_multiqc": "stat",
          "star_multiqc": "logs",
          "bam_star_rsem": "bam_star",
          "bam_genome_rsem": "bam_genome",
          "bam_transcript_rsem": "bam_transcript"
        ]
        directives:
          - label: [ midmem, midcpu ]
      )
    
      // RSEM_Star BAM
      | samtools_sort.run(
        runIf: { id, state -> state.aligner == 'star_rsem' },
        fromState: ["input": "bam_star_rsem"],
        toState: ["genome_bam_sorted": "output"],
        key: "genome_sorted"
        directives:
          - label: [ midmem, midcpu ]
      )
      | samtools_index.run(
        runIf: { id, state -> state.aligner == 'star_rsem' },
        fromState: [ 
          "input": "genome_bam_sorted", 
          "csi": "bam_csi_index"
        ],
        toState: [ "genome_bam_index": "output" ],
        key: "genome_sorted"
        directives:
          - label: [ midmem, midcpu ]
      )
      | samtools_stats.run(
        runIf: { id, state -> state.aligner == 'star_rsem' },
        fromState: [
          "input": "genome_bam_sorted", 
          "bai": "genome_bam_index", 
          "fasta": "fasta"
        ],
        toState: [ "genome_bam_stats": "output" ],
        key: "genome_stats"
        directives:
          - label: [ midmem, midcpu ]
      )
      | samtools_flagstat.run(
        runIf: { id, state -> state.aligner == 'star_rsem' },
        fromState: [
          "bam": "genome_bam_sorted", 
          "bai": "genome_bam_index", 
          "fasta": "fasta" 
        ],
        toState: [ "genome_bam_flagstat": "output" ],
        key: "genome_flagstat"
        directives:
          - label: [ midmem, midcpu ]
      )
      | samtools_idxstats.run(
        runIf: { id, state -> state.aligner == 'star_rsem' },
        fromState: [
          "bam": "genome_bam_sorted", 
          "bai": "genome_bam_index", 
          "fasta": "fasta" 
        ],
        toState: [ "genome_bam_idxstats": "output" ],
        key: "genome_idxstats"
        directives:
          - label: [ midmem, midcpu ]
      )
      | map { id, state -> 
        def mod_state = state.findAll { key, value -> value instanceof java.nio.file.Path && value.exists() }
        [ id, mod_state ]
      }

      | setState (
        [ "star_multiqc": "star_multiqc", 
          "rsem_multiqc": "rsem_multiqc",
          "salmon_multiqc": "salmon_multiqc",
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
          "quant_out_dir": "quant_out_dir",
          "quant_results_file": "quant_results_file",
          "rsem_counts_gene": "rsem_counts_gene",
          "rsem_counts_transcripts": "rsem_counts_transcripts",
          "bam_genome_rsem": "bam_genome_rsem",
          "bam_transcript_rsem": "bam_transcript_rsem" ]
      )

  emit:
    output_ch
}
