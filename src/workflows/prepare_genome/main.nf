workflow run_wf {

    take: 
        input_ch

    main: 
        output_ch = input_ch 

        // decompress fasta
        | gunzip.run (
            fromState: ["input": "fasta"], 
            toState: ["fasta": "output"], 
            key: "gunzip_fasta" 
        )

        // decompress gtf or gff
        | map { id, state ->
            def input = state.gtf ? state.gtf : state.gff
            [ id, state + [input: input] ]
        }
        | gunzip.run ( 
            fromState: ["input": "input"], 
            toState: ["annotation": "output"], 
            key: "gunzip_annotation"
        )

        | map { id, state ->
            def gtf = state.gtf ? state.annotation : ''
            [ id, state + [gtf: gtf] ]
        }

        // gff to gtf
        | gffread.run (
            runIf: {id, state -> state.gff}, 
            fromState: ["input": "annotation"], 
            toState: ["gtf": "output"] 
        )

        // decompress additional fasta
        | gunzip.run (
            runIf: {id, state -> state.additional_fasta}, 
            fromState: ["input": "additional_fasta"], 
            toState: ["additional_fasta": "output"], key: "gunzip_additional_fasta" 
        )

        // concatenate additional fasta
        | cat_additional_fasta.run (
            runIf: {id, state -> state.additional_fasta}, 
            fromState: [
                "fasta": "fasta", 
                "gtf": "gtf", 
                "additional_fasta": "additional_fasta", 
                "biotype": "biotype" ], 
            toState: [
                "fasta": "fasta_output", 
                "gtf": "gtf_output" ], 
            key: "cat_additional" 
        ) 

        // decompress bed file
        | gunzip.run (
            runIf: {id, state -> state.gene_bed}, 
            fromState: ["input": "gene_bed"], 
            toState: ["gene_bed": "output"], 
            key: "gunzip_gene_bed" 
        )

        // gtf to bed 
        | gtf2bed.run (
            runIf: { id, state -> !state.gene_bed}, 
            fromState: ["gtf": "gtf"], 
            toState: ["gene_bed": "bed_output"] 
        ) 

        // decompress transcript fasta
        | gunzip.run (
            runIf: {id, state -> state.transcript_fasta}, 
            fromState: ["input": "transcript_fasta"], 
            toState: ["transcript_fasta": "output"], 
            key: "transcript_fasta"
        )

        // preprocess transcripts fasta if gtf is in gencode format
        | preprocess_transcripts_fasta.run (
            runIf: {id, state -> state.transcript_fastaa && state.gencode}, 
            fromState: ["transcript_fasta": "transcript_fasta"], 
            toState: ["transcript_fasta": "output"], 
            key: "transcript_fixed" 
        )

        // filter gtf for genes in genome
        | gtf_gene_filter.run (
            runIf: {id, state -> !state.transcript_fasta}, 
            fromState: ["fasta": "fasta", "gtf": "gtf"], 
            toState: ["filtered_gtf": "filtered_gtf"] 
        )

        // prepare reference for RSEM
        | rsem_prepare_reference.run ( 
            runIf: {id, state -> !state.transcript_fasta}, 
            fromState: [
                "fasta": "fasta", 
                "gtf": "filtered_gtf" ], 
            toState: ["transcript_fasta": "transcript_fasta"], 
            key: "rsem_ref" 
        )

        // chromosome size and fai index
        | getchromsize.run (
            fromState: ["fasta": "fasta"], 
            toState: [
                "fai": "fai", 
                "sizes": "sizes" ], 
            key: "chromsizes" 
        )
        
        // untar bbsplit index, if available
        | untar.run (
            runIf: {id, state -> state.bbsplit_index}, 
            fromState: ["input": "bbsplit_index"], 
            toState: ["bbsplit_index": "output"], 
            key: "bbsplit_uncompressed" 
        )
        
        // create bbsplit index, if not already availble
        | bbmap_bbsplit.run (
            runIf: {id, state -> !state.bbsplit_index}, 
            fromState: [ 
                "primary_ref": "fasta", 
                "bbsplit_fasta_list": "bbsplit_fasta_list" ], 
            toState: ["bbsplit_index": "bbsplit_index"], 
            args: ["only_build_index": true], 
            key: "bbsplit_index_uncompressed" 
        )

        // Uncompress STAR index or generate from scratch if required
        | untar.run (
            runIf: {id, state -> state.star_index}, 
            fromState: ["input": "star_index"], 
            toState: ["star_index": "output"], 
            key: "star_index_uncompressed"
        )
        
        | star_genomegenerate.run (
            runIf: {id, state -> !state.star_index}, 
            fromState: [ 
                "fasta": "fasta", 
                "gtf": "gtf" ], 
            toState: ["star_index": "star_index"], 
            key: "star_index_uncompressed"
        )

        // TODO: Uncompress RSEM index or generate from scratch if required

        // TODO: Uncompress HISAT2 index or generate from scratch if required

        // Uncompress Salmon index or generate from scratch if required
        | untar.run (
            runIf: {id, state -> state.salmon_index}, 
            fromState: ["input": "salmon_index"], 
            toState: ["salmon_index": "output"], 
            key: "salmon_index_uncompressed"
        )

        | salmon_index.run (
            runIf: {id, state -> !state.salmon_index}, 
            fromState: [ 
                "genome_fasta": "fasta", 
                "transcriptome_fasta": "transcript_fasta" ], 
            toState: ["salmon_index": "salmon_index"], 
            key: "salmon_index_uncompressed" 
        )

        | setState ( 
            [ "fasta_uncompressed": "fasta", 
            "gtf_uncompressed": "gtf", 
            "transcript_fasta_uncompressed": "transcript_fasta", 
            "star_index_uncompressed": "star_index", 
            "salmon_index_uncompressed": "salmon_index", 
            "bbsplit_index_uncompressed": "bbsplit_index", 
            "chrom_sizes": "sizes", 
            "fai": "fai" ] )      

    emit: 
        output_ch
}
