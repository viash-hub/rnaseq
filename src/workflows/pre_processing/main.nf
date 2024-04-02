// TODO: Add arguments for fastp

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
        input: input, 
        versions: state.versions ]
      },
      toState: [
        "fastqc_html_1": "fastqc_html_1",
        "fastqc_html_2": "fastqc_html_2",
        "fastqc_zip_1": "fastqc_zip_1",
        "fastqc_zip_2": "fastqc_zip_2", 
        "versions": "updated_versions" 
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
        umi_discard_read: state.umi_discard_read, 
        versions: state.versions ]
      },
      toState: [ 
        "fastq_1": "fastq_1", 
        "fastq_2": "fastq_2", 
        "versions": "updated_versions" 
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
        extra_trimgalore_args: state.extra_trimgalore_args, 
        min_trimmed_reads: state.min_trimmed_reads, 
        versions: state.versions ]
      },
      toState: [
        "fastq_1": "fastq_1", 
        "fastq_2": "fastq_2",
        "trim_log_1": "trim_log_1", 
        "trim_log_2": "trim_log_2", 
        "trim_zip_1": "trim_zip_1",
        "trim_zip_2": "trim_zip_2",
        "trim_html_1": "trim_html_1",
        "trim_html_2": "trim_html_2", 
        "versions": "updated_versions" 
      ]
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
        input: input, 
        versions: state.versions ]
      },
      toState: [
        "trim_html_1": "fastqc_html_1",
        "trim_html_2": "fastqc_html_2",
        "trim_zip_1": "fastqc_zip_1",
        "trim_zip_2": "fastqc_zip_2", 
        "versions": "updated_versions" 
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
        built_bbsplit_index: state.bbsplit_index,
        versions: state.versions ]
      },
      args: ["only_build_index": false], 
      toState: [
        "fastq_1": "fastq_1", 
        "fastq_2": "fastq_2",
        "versions": "updated_versions"
      ]
    )

    // Sort reads by rRNA and non-rRNA
    | sortmerna.run (
      runIf: { id, state -> state.remove_ribo_rna },
      fromState: { id, state ->
        def input = state.paired ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
        [ paired: state.paired,
          input: input,
          ribo_database_manifest: state.ribo_database_manifest, 
          versions: state.versions ] 
      },
      toState: [
        "fastq_1": "fastq_1", 
        "fastq_2": "fastq_2",
        "sortmerna_log": "sortmerna_log", 
        "versions": "updated_versions" 
      ] 
    )

    // Sub-sample FastQ files and pseudo-align with Salmon to auto-infer strandedness
    | fq_subsample.run (
      runIf: { id, state -> state.strandedness == 'auto' }, 
      fromState: { id, state ->
        def input = state.paired ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
        [ input: input,
          extra_args: state.extra_fq_subsample_args, 
          versions: state.versions ] 
      },
      toState: [
        "subsampled_fastq_1": "output_1",
        "subsampled_fastq_2": "output_2", 
        "versions": "updated_versions" 
      ]
    )

    | salmon_quant.run (
      runIf: { id, state -> state.strandedness == 'auto' }, 
      fromState: { id, state ->
        def input = state.paired ? [ state.subsampled_fastq_1, state.subsampled_fastq_2 ] : [ state.subsampled_fastq_1 ]
        [ paired: state.paired, 
          strandedness: state.strandedness, 
          input: input, 
          transcript_fasta: state.transcript_fasta, 
          gtf: state.gtf, 
          salmon_index: state.salmon_index, 
          versions: state.versions ]
      },
      args: [
        "alignment_mode": false, 
        "lib_type": "A", 
        "extra_salmon_quant_args": "--skipQuant"
      ],
      toState: [
        "salmon_quant_output": "output",
        "salmon_json_info": "json_info", 
        "versions": "updated_versions" 
      ]
    )

    | map { id, state -> 
      def mod_state = (!state.paired) ? 
        [trim_log_2: state.remove(state.trim_log_2), trim_zip_2: state.remove(state.trim_zip_2), trim_html_2: state.remove(state.trim_html_2), failed_trim_unpaired2: state.remove(state.failed_trim_unpaired2)] : 
        []
      [ id, state + mod_state ]
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
        "salmon_json_info": "salmon_json_info", 
        "failed_trim": "failed_trim",
        "failed_trim_unpaired1": "failed_trim_unpaired1",
        "failed_trim_unpaired2": "failed_trim_unpaired2",
        "trim_json": "trim_json",
        "trim_html": "trim_html",
        "trim_merged_out": "trim_merged_out",
        "updated_versions": "versions" 
    )

  emit:
    output_ch
}
