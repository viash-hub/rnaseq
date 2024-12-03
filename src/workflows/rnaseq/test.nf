include { rnaseq } from params.rootDir + "/target/nextflow/workflows/rnaseq/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

  resources_test = file(params.resources_test)

  output_ch = 

    Channel.fromList([
      [
        id: "_test",
        fastq_1: "https://github.com/nf-core/test-datasets/raw/refs/heads/rnaseq/testdata/GSE110004/SRR6357070_1.fastq.gz;" +
          "https://github.com/nf-core/test-datasets/raw/refs/heads/rnaseq/testdata/GSE110004/SRR6357071_1.fastq.gz" ,
        fastq_2: "https://github.com/nf-core/test-datasets/raw/refs/heads/rnaseq/testdata/GSE110004/SRR6357070_2.fastq.gz;" +
          "https://github.com/nf-core/test-datasets/raw/refs/heads/rnaseq/testdata/GSE110004/SRR6357071_2.fastq.gz",
        strandedness: "reverse",
        fasta: "https://github.com/nf-core/test-datasets/raw/refs/heads/rnaseq/reference/genome.fasta",
        gtf: "https://github.com/nf-core/test-datasets/raw/refs/heads/rnaseq/reference/genes.gtf.gz",
        additional_fasta: "https://github.com/nf-core/test-datasets/raw/refs/heads/rnaseq/reference/gfp.fa.gz",
        transcript_fasta: "https://github.com/nf-core/test-datasets/raw/refs/heads/rnaseq/reference/transcriptome.fasta",
        salmon_index: "https://github.com/nf-core/test-datasets/raw/refs/heads/rnaseq/reference/salmon.tar.gz",
        publish_dir: "output/",
        skip_deseq2_qc: true
      ]
    ])

    | map { state -> [state.id, state] }

    | rnaseq.run(
      fromState: { id, state -> state },
      toState: { id, output, state -> output }
    )

    | view { output ->
      assert output.size() == 2 : "Outputs should contain two elements; [id, state]"

      // check id
      def id = output[0]
      assert id.endsWith("_test")

      // check output
      def state = output[1]
      assert state instanceof Map : "State should be a map. Found: ${state}"
      assert state.containsKey("output_fasta") : "Output should contain key 'output_fasta'."
      assert state.output_fasta.isFile() : "'output_fasta' should be a file."
      assert state.containsKey("output_gtf") : "Output should contain key 'output_gtf'."
      assert state.output_gtf.isFile() : "'output_gtf' should be a file."

      "Output: $output"
    }
}
