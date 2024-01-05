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

    // Sub-sample FastQ files and pseudo-align with Salmon to auto-infer strandedness
    | fq_subsample.run (
      runIf: { id, state -> state.strandedness == 'auto' }, 
      fromState: [ 
        "input": "input", 
        "extra_args": "extra_fq_subsample_args", 
        "versions": "versions" 
      ],
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

    // Trim reads
    | trimgalore.run (
      runIf: { id, state -> !state.skip_trimming },
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

    // Filter out contaminant RNA
    | bbmap_bbsplit.run (
      runIf: { id, state -> !state.skip_bbsplit },
      fromState: { id, state ->
        def input = state.paired ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
        [ paired: state.paired,
        input: input,
        built_bbsplit_index: state.bbsplit_index,
        bbsplit_fasta_list: state.bbsplit_fasta_list,
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

    | map { id, state -> 
      def mod_state = (!state.paired) ? 
        [trim_log_2: state.remove(state.trim_log_2), trim_zip_2: state.remove(state.trim_zip_2), trim_html_2: state.remove(state.trim_html_2)] : 
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
        "updated_versions": "versions" 
    )

  emit:
    output_ch
}
