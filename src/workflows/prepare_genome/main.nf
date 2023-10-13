workflow run_wf {

    take: 
        input_ch

    main: 
        output_ch = input_ch 
        | view {"Input: $it"}

        // decompress fasta
        | gunzip.run(auto: [publish: true], fromState: ["input": "fasta"], toState: ["fasta_uncompressed": "output"], key: "gunzip_fasta")
        | view {"State: $it"}

        // decompress gtf or gff
        | map { id, state ->
        def input = state.gtf ? state.gtf : state.gff
        [ id, state + [input: input] ]
        }
        | gunzip.run(auto: [publish: true], fromState: ["input": "input"], toState: ["annotation": "output"], key: "gunzip_annotation")
        | view {"State: $it"}

        | map { id, state ->
        def gtf = state.gtf ? state.annotation : ''
        [ id, state + [gtf_uncompressed: gtf] ]
        }

        // gff to gtf
        | gffread.run(runIf: {id, state -> state.gff}, auto: [publish: true], fromState: ["input": "annotation"], toState: ["gtf_uncompressed": "output"])

        // decompress additional fasta
        | gunzip.run(runIf: {id, state -> state.additional_fasta}, auto: [publish: true], fromState: ["input": "additional_fasta"], toState: ["additional_fasta": "output"], key: "gunzip_additional_fasta")

        // concatenate additional fasta
        | cat_additional_fasta.run(runIf: {id, state -> state.additional_fasta}, auto: [publish: true], fromState: ["fasta": "fasta_uncompressed", "gtf": "gtf_uncompressed", "additional_fasta": "additional_fasta", "biotype": "biotype"], toState: ["fasta_uncompressed": "fasta_output", "gtf_uncompressed": "gtf_output"], key: "cat_additional") 

        // decompress bed file
        | gunzip.run(runIf: {id, state -> state.gene_bed}, auto: [publish: true], fromState: ["input": "gene_bed"], toState: ["gene_bed_uncompressed": "output"],key: "gunzip_gene_bed")

        // gtf to bed (not working?)
        // | gtf2bed.run(runIf: { id, state -> !state.gene_bed}, auto: [publish: true], fromState: ["gtf": "gtf_uncompressed"], toState: ["gene_bed_uncompressed": "bed_output"]) 
        // | view {"State: $it"}

        // decompress transcript fasta
        | gunzip.run(runIf: {id, state -> state.transcript_fast}, auto: [publish: true], fromState: ["input": "transcript_fasta"], toState: ["transcript_fasta_uncompressed": "output"], key: "transcript_fasta")

        // preprocess transcripts fasta if gtf is in gencode format
        | preprocess_transcripts_fasta.run(runIf: {id, state -> state.transcript_fastaa && state.gencode}, auto: [publish: true], fromState: ["transcript_fasta": "transcript_fasta_uncompressed"], toState: ["transcript_fasta_uncompressed": "output"], key: "transcript_fixed")

        // filter gtf for genes in genome
        | gtf_gene_filter.run(runIf: {id, state -> !state.transcript_fasta}, auto: [publish: true], fromState: ["fasta": "fasta_uncompressed", "gtf": "gtf_uncompressed"], toState: ["filtered_gtf": "filtered_gtf"])

        // prepare reference for RSEM
        | rsem_prepare_reference.run(runIf: {id, state -> !state.transcript_fasta}, auto: [publish: true], fromState: ["fasta": "fasta_uncompressed", "gtf": "filtered_gtf"], toState: ["transcript_fasta_uncompressed": "transcript_fasta"], key: "rsem_ref")

        // chromosome size and fai index
        | getchromsize.run(auto: [publish: true], fromState: ["fasta": "fasta_uncompressed"], toState: ["fai": "fai", "sizes": "sizes"], key: "chromsizes")
        
        // untar bbsplit index, if available
        | untar.run(runIf: {id, state -> state.bbsplit_index}, auto: [publish: true], fromState: ["input": "bbsplit_index"], toState: ["bbsplit_index_uncompressed": "output"], key: "bbsplit_uncompressed")
        
        // create bbsplit index, if not already availble
        // | map { id, state -> [ id, state + [ "only_build_index": true ] ] }
        | bbmap_bbsplit.run(runIf: {id, state -> !state.bbsplit_index}, auto: [publish: true], fromState: ["primary_ref": "fasta_uncompressed", "bbsplit_fasta_list": "bbsplit_fasta_list"], toState: ["bbsplit_index_uncompressed": "bbsplit_index"], args: ["only_build_index": true], key: "bbsplit_index_uncompressed")

        // Uncompress STAR index or generate from scratch if required
        | untar.run(runIf: {id, state -> state.star_index}, auto: [publish: true], fromState: ["input": "star_index"], toState: ["star_index_uncompressed": "output"], key: "star_uncompressed")
        
        | star_genomegenerate.run(runIf: {id, state -> !state.star_index}, auto: [publish: true], fromState: ["fasta": "fasta_uncompressed", "gtf": "gtf_uncompressed"], toState: ["star_index_uncompressed": "star_index"], key: "star_uncompressed")

        // TODO: Uncompress RSEM index or generate from scratch if required

        // TODO: Uncompress HISAT2 index or generate from scratch if required

        // Uncompress Salmon index or generate from scratch if required
        | untar.run(runIf: {id, state -> state.salmon_index}, auto: [publish: true], fromState: ["input": "salmon_index"], toState: ["salmon_index_uncompressed": "output"], key: "salmon_index_uncompressed")

        | salmon_index.run(runIf: {id, state -> !state.salmon_index}, auto: [publish: true], fromState: ["genome_fasta": "fasta_uncompressed", "transcriptome_fasta": "transcript_fasta_uncompressed"], toState: ["salmon_index_uncompressed": "salmon_index"], key: "salmon_index_uncompressed")

        | view {"Output: $it"}  

      //   | setState([
      //   "fasta_uncompressed": "fasta", 
      //   "gtf_uncompressed": "gtf",
      //   "transcript_fasta_uncompressed": "transcript_fasta", 
      //   "star_index_uncompressed": "star_index" 
      //   "salmon_index_uncompressed": "salmon_index", 
      //   "bbsplit_index_uncompressed": "bbsplit_index", 
      //   "chrom_sizes": "sizes", 
      //   "fai": "fai"
      // ])      

    emit: 
        output_ch
}
