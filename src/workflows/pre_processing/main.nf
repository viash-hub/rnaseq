workflow run_wf {
  
  take:
    input_ch

  main:
    output_ch = input_ch
    
    | map { id, state ->
      def input = state.fastq_2 ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
      def paired = input.size() == 2
      [ id, state + [paired: paired, input: input] ]
    }

    // Perform QC on input fastq files
    | fastqc.run (
      runIf: { id, state -> !state.skip_qc && !state.skip_fastqc },
      fromState: { id, state ->
        def input = state.paired ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
        [ paired: state.paired,
        input: input ]
      },
      toState: [
        "fastqc_html_1": "fastqc_html_1",
        "fastqc_html_2": "fastqc_html_2",
        "fastqc_zip_1": "fastqc_zip_1",
        "fastqc_zip_2": "fastqc_zip_2"
      ]
    )

    // Extract UMIs from fastq files and discard read 1 or read 2 if required
    | umitools_extract.run (
      runIf: { id, state -> state.with_umi && !state.skip_umi_extract },
      fromState: { id, state ->
        def input = state.paired ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
        def bc_pattern = state.paired ? [ state.umitools_bc_pattern, state.umitools_bc_pattern2 ] : [ state.umitools_bc_pattern ]
        [ paired: state.paired,
        input: input, 
        bc_pattern: bc_pattern, 
        umi_discard_read: state.umi_discard_read ]
      },
      toState: [ 
        "fastq_1": "fastq_1", 
        "fastq_2": "fastq_2"
      ]
    )
    
    // Discard read if required
    | map { id, state -> 
      def paired = state.paired
      def fastq_2 = state.fastq_2
      if (paired && state.with_umi && !state.skip_umi_extract && state.umi_discard_read != 0) {
        fastq_2 = state.remove(state.fastq_2) 
        paired = false
      }
      [ id, state + [paired: paired, fastq_2: fastq_2] ]
    }

    // Trim reads using Trim galore!
    | trimgalore.run (
      runIf: { id, state -> !state.skip_trimming && state.trimmer == "trimgalore" },
      fromState: { id, state ->
        def input = state.paired ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
        [ paired: state.paired,
        input: input, 
        min_trimmed_reads: state.min_trimmed_reads ]
      },
      toState: [
        "fastq_1": "trimmed_r1", 
        "fastq_2": "trimmed_r2",
        "trim_log_1": "trimming_report_r1", 
        "trim_log_2": "trimming_report_r2", 
        "trim_zip_1": "trimmed_fastqc_zip_1",
        "trim_zip_2": "trimmed_fastqc_zip_2",
        "trim_html_1": "trimmed_fastqc_html_1",
        "trim_html_2": "trimmed_fastqc_html_2"
      ],
      args: [gzip: true]
    )

    // Trim reads using fastp
    | fastp.run(
      runIf: { id, state -> !state.skip_trimming && state.trimmer == "fastp" },
      fromState: [
        "in1": "fastq_1",
        "in2": "fastq_2",
        "merge": "fastp_save_merged", 
        "interleaved_in": "interleaved_reads",
        "detect_adapter_for_pe": "fastp_pe_detect_adapter",
        "adapter_fasta": "fastp_adapter_fasta"
      ],
      toState: [
        "fastq_1": "out1", 
        "fastq_2": "out2",
        "failed_trim": "failed_out",
        "failed_trim_unpaired1": "unpaired1",
        "failed_trim_unpaired2": "unpaired2",
        "trim_json": "json",
        "trim_html": "html",
        "trim_merged_out": "merged_out"
      ]
    )

    // Perform FASTQC on reads trimmed using fastp
    | fastqc.run(
      runIf: { id, state -> !state.skip_trimming && state.trimmer == "fastp" },
      fromState: { id, state ->
        def input = state.paired ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
        [ paired: state.paired,
        input: input ]
      },
      toState: [
        "trim_html_1": "fastqc_html_1",
        "trim_html_2": "fastqc_html_2",
        "trim_zip_1": "fastqc_zip_1",
        "trim_zip_2": "fastqc_zip_2"
      ], 
      key: "fastqc_trimming"
    )

    // Filter out contaminant RNA
    | bbmap_bbsplit.run (
      runIf: { id, state -> !state.skip_bbsplit },
      fromState: { id, state ->
        def input = state.paired ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
        [ paired: state.paired,
        input: input,
        built_bbsplit_index: state.bbsplit_index ]
      },
      args: ["only_build_index": false], 
      toState: [
        "fastq_1": "fastq_1", 
        "fastq_2": "fastq_2"
      ]
    )

    // Sort reads by rRNA and non-rRNA
    | sortmerna.run (
      runIf: { id, state -> state.remove_ribo_rna },
      fromState: { id, state ->
        def input = state.paired ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
        [ paired: state.paired,
          input: input,
          ribo_database_manifest: state.ribo_database_manifest ] 
      },
      toState: [
        "fastq_1": "fastq_1", 
        "fastq_2": "fastq_2",
        "sortmerna_log": "sortmerna_log"
      ] 
    )

    // Sub-sample FastQ files and pseudo-align with Salmon to auto-infer strandedness
    | fq_subsample.run (
      runIf: { id, state -> state.strandedness == 'auto' }, 
      fromState: { id, state ->
        def input = state.paired ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
        [ 
          input: input,
          extra_args: state.extra_fq_subsample_args
        ]
      },
      toState: [
        "subsampled_fastq_1": "output_1",
        "subsampled_fastq_2": "output_2"
      ]
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
    | salmon_quant.run (
      runIf: { id, state -> state.strandedness == 'auto' }, 
      fromState: { id, state ->
        def unmated_reads = !state.paired ? state.subsampled_fastq_1 : null
        def mates1 = state.paired ? state.subsampled_fastq_1 : null
        def mates2 = state.paired ? state.subsampled_fastq_2 : null
        [ unmated_reads: unmated_reads,
          mates1: mates1, 
          mates2: mates2, 
          gene_map: state.gtf, 
          index: state.salmon_index,
          lib_type: state.lib_type ]
      },
      args: [ "skip_quant": true ],
      toState: [ "salmon_quant_output": "output" ]
    )

    | map { id, state -> 
      def mod_state = (!state.paired) ? 
        [trim_log_2: state.remove(state.trim_log_2), trim_zip_2: state.remove(state.trim_zip_2), trim_html_2: state.remove(state.trim_html_2), failed_trim_unpaired2: state.remove(state.failed_trim_unpaired2)] : 
        []
      [ id, state + mod_state ]
    }

    | map { id, state -> 
      def mod_state = state.findAll { key, value -> value instanceof java.nio.file.Path && value.exists() }
      [ id, mod_state ]
    }

    | setState ( 
        "fastqc_html_1": "fastqc_html_1",
        "fastqc_html_2": "fastqc_html_2",
        "fastqc_zip_1": "fastqc_zip_1",
        "fastqc_zip_2": "fastqc_zip_2", 
        "qc_output1": "fastq_1",
        "qc_output2": "fastq_2", 
        "trim_log_1": "trim_log_1", 
        "trim_log_2": "trim_log_2", 
        "trim_zip_1": "trim_zip_1",
        "trim_zip_2": "trim_zip_2",
        "trim_html_1": "trim_html_1",
        "trim_html_2": "trim_html_2",
        "sortmerna_log": "sortmerna_log",
        "failed_trim": "failed_trim",
        "failed_trim_unpaired1": "failed_trim_unpaired1",
        "failed_trim_unpaired2": "failed_trim_unpaired2",
        "trim_json": "trim_json",
        "trim_html": "trim_html",
        "trim_merged_out": "trim_merged_out",
        "salmon_quant_output": "salmon_quant_output"
    )

  emit:
    output_ch
}
