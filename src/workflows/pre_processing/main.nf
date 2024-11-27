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

    | fastqc.run (
      runIf: { id, state -> !state.skip_qc && !state.skip_fastqc },
      fromState: [ "input": "input" ], 
      toState: {id, output_state, state -> 
          def newKeys = [
            "fastqc_html_1":output_state["html"][0],
            "fastqc_html_2": output_state["html"][1],
            "fastqc_zip_1": output_state["zip"][0],
            "fastqc_zip_2": output_state["zip"][1]
          ]
          def new_state = state + newKeys
          return new_state
        },
      args: [html: "*.html", zip: "*.zip"]
    )

    // Extract UMIs from fastq files and discard read 1 or read 2 if required
    | umi_tools_extract.run (
      runIf: { id, state -> state.with_umi && !state.skip_umi_extract },
      fromState: { id, state ->
        def bc_pattern2 = state.paired ? state.umitools_bc_pattern2 : state.remove(state.umitools_bc_pattern2)
        def output = "${id}.r1.fastq.gz"
        def read2_out = state.paired ? "${id}.r2.fastq.gz" : state.remove(state.fastq_2)
        [ input: state.fastq_1, 
        read2_in: state.fastq_2,
        bc_pattern: state.umitools_bc_pattern, 
        bc_pattern2: bc_pattern2,
        extract_method: state.umitools_extract_method, 
        umi_separator: state.umitools_umi_separator,
        grouping_method: state.umitools_grouping_method,
        output: output,
        read2_out: read2_out ]
      },
      toState: [ 
        "fastq_1": "output", 
        "fastq_2": "read2_out"
      ]
    )
    
    // Discard read if required
    | map { id, state -> 
      def paired = state.paired
      def fastq_1 = state.fastq_1
      def fastq_2 = state.fastq_2
      if (paired && state.with_umi && !state.skip_umi_extract && state.umi_discard_read != 0) {
        if (state.umi_discard_read == 1) {
          fastq_1 = fastq_2
        } 
        fastq_2 = state.remove(state.fastq_2)
        paired = false
      }
      [ id, state + [paired: paired, fastq_1: fastq_1, fastq_2: fastq_2] ]
    }

    // Trim reads using Trim galore!
    | trimgalore.run (
      runIf: { id, state -> !state.skip_trimming && state.trimmer == "trimgalore" },
      fromState: { id, state ->
        def input = state.paired ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
        [ paired: state.paired,
        input: input, 
        min_trimmed_reads: state.min_trimmed_reads, 
        trimmed_r1: state.qc_output1, 
        trimmed_r2: state.qc_output2 ]
      },
      args: [gzip: true, fastqc: true],
      toState: [
        "fastq_1": "trimmed_r1", 
        "fastq_2": "trimmed_r2",
        "trim_log_1": "trimming_report_r1", 
        "trim_log_2": "trimming_report_r2", 
        "trim_zip_1": "trimmed_fastqc_zip_1",
        "trim_zip_2": "trimmed_fastqc_zip_2",
        "trim_html_1": "trimmed_fastqc_html_1",
        "trim_html_2": "trimmed_fastqc_html_2"
      ]
    )

    // Trim reads using fastp
    | fastp.run(
      runIf: { id, state -> !state.skip_trimming && state.trimmer == "fastp" },
      fromState: { id, state -> 
        def outputState = state.paired ? [out1: state.qc_output1, out2: state.qc_output2] : [out1: state.qc_output1, out2: state.remove(state.qc_output2)]
        [input_1: state.fastq_1, input_2: state.fastq_2] + outputState
        [ in1: state.fastq_1,
        in2: state.fastq_2,
        merge: state.fastp_save_merged, 
        interleaved_in: state.interleaved_reads,
        detect_adapter_for_pe: state.paired,
        adapter_fasta: state.fastp_adapter_fasta ] + outputState
      },
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
    | fastqc.run (
      runIf: { id, state -> !state.skip_trimming && state.trimmer == "fastp" },
      fromState: { id, state ->
        def input = state.paired ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
        [ input: input ]
      }, 
      toState: {id, output_state, state -> 
          def newKeys = [
            "trim_html_1":output_state["html"][0],
            "trim_html_2": output_state["html"][1],
            "trim_zip_1": output_state["zip"][0],
            "trim_zip_2": output_state["zip"][1]
          ]
          def new_state = state + newKeys
          return new_state
        },
      args: [html: "*.html", zip: "*.zip"],
      key: "fastqc_trimming"
    )

    // Filter out contaminant RNA
    | bbmap_bbsplit.run (
      runIf: { id, state -> !state.skip_bbsplit },
      fromState: { id, state ->
        def input = state.paired ? [ state.fastq_1, state.fastq_2 ] : [ state.fastq_1 ]
        [ paired: state.paired,
        input: input,
        build: state.bbsplit_index ]
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
        def filePaths = state.ribo_database_manifest.readLines()
        def refs = filePaths.collect { it }
        def other = "${id}_non_rRNA_reads/"
        [ paired_in: state.paired,
          input: input,
          ref: refs,
          out2: state.paired, 
          other: other ] 
      },
      args: [fastx: true, num_alignments: 1],
      toState: { id, output_state, state -> 
        def newKeys = [ 
          "sortmerna_output": output_state["other"],
          "sortmerna_log": output_state["log"]
        ]
        def new_state = state + newKeys
        return new_state
      }
    )
    | map { id, state -> 
      if (state.remove_ribo_rna) {
        def fastq_1 = state.sortmerna_output.listFiles().find{it.name == "other_fwd.fq.gz"}
        def fastq_2 = state.sortmerna_output.listFiles().find{it.name == "other_rev.fq.gz"}
        [ id, state + [fastq_1: fastq_1, fastq_2: fastq_2] ] 
      } else {
        [ id, state ]
      }   
    }

    // Sub-sample FastQ files and pseudo-align with Salmon to auto-infer strandedness
    | fq_subsample.run (
      runIf: { id, state -> state.strandedness == 'auto' }, 
      fromState: { id, state -> 
        def outputState = state.paired ? [output_1: state.qc_output1, output_2: state.qc_output2] : [output_1: state.qc_output1, output_2: state.remove(state.qc_output2)]
        [input_1: state.fastq_1, input_2: state.fastq_2] + outputState
      },
      args: [
        record_count: 1000, 
        seed: 1
      ],
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
        "trim_json": "trim_json",
        "trim_html": "trim_html",
        "trim_merged_out": "trim_merged_out",
        "salmon_quant_output": "salmon_quant_output"
    )

  emit:
    output_ch
}
