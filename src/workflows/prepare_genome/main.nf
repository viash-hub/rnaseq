workflow run_wf {

    take: 
        input_ch

    main: 
        output_ch = input_ch 

        // decompress fasta
        | gunzip.run (
            fromState: ["input": "fasta"], 
            toState: ["fasta": "output"], 
            key: "gunzip_fasta",
            args: [output: "reference_genome.fasta"] 
        )

        // decompress gtf
        | gunzip.run ( 
            runIf: {id, state -> state.gtf},
            fromState: ["input": "gtf"], 
            toState: ["gtf": "output"], 
            key: "gunzip_gtf",
            args: [output: "gene_annotation.gtf"]
        )

        // decompress gff
        | gunzip.run ( 
            runIf: {id, state -> !state.gtf && state.gff},
            fromState: ["input": "gff"], 
            toState: ["gff": "output"], 
            key: "gunzip_gff",
            args: [output: "gene_annotation.gff"] 
        )

        // gff to gtf
        | gffread.run (
            runIf: {id, state -> state.gff}, 
            fromState: ["input": "annotation"], 
            toState: ["gtf": "output"],
            args: [output: "gene_annotation.gtf"] 
        )

        // decompress additional fasta
        | gunzip.run (
            runIf: {id, state -> state.additional_fasta}, 
            fromState: ["input": "additional_fasta"], 
            toState: ["additional_fasta": "output"], 
            key: "gunzip_additional_fasta",
            args: [output: "additional.fasta"]  
        )

        // concatenate additional fasta
        | cat_additional_fasta.run (
            runIf: {id, state -> state.additional_fasta}, 
            fromState: [
                "fasta": "fasta", 
                "gtf": "gtf", 
                "additional_fasta": "additional_fasta", 
                "biotype": "biotype" 
            ], 
            toState: [
                "fasta": "fasta_output", 
                "gtf": "gtf_output" 
            ], 
            key: "cat_additional",
            args: [
                fasta_output: "genome_additional_concat.fasta", 
                gtf_output: "genome_additional_concat.gtf"]  
        ) 

        // decompress bed file
        | gunzip.run (
            runIf: {id, state -> state.gene_bed}, 
            fromState: ["input": "gene_bed"], 
            toState: ["gene_bed": "output"], 
            key: "gunzip_gene_bed",
            args: [output: "genome_additional_concat.bed"]
        )

        // gtf to bed 
        | gtf2bed.run (
            runIf: { id, state -> !state.gene_bed}, 
            fromState: ["gtf": "gtf"], 
            toState: ["gene_bed": "bed_output"], 
            args: [bed_output: "genome_additional_concat.bed"]
        ) 

        // decompress transcript fasta
        | gunzip.run (
            runIf: {id, state -> state.transcript_fasta}, 
            fromState: ["input": "transcript_fasta"], 
            toState: ["transcript_fasta": "output"], 
            key: "transcript_fasta", 
            args: [output: "transcriptome.fasta"]
        )

        // preprocess transcripts fasta if gtf is in gencode format
        | preprocess_transcripts_fasta.run (
            runIf: {id, state -> state.transcript_fastaa && state.gencode}, 
            fromState: ["transcript_fasta": "transcript_fasta"], 
            toState: ["transcript_fasta": "output"], 
            key: "transcript_fixed", 
            args: [output: "transcriptome.fasta"] 
        )

        // filter gtf for genes in genome
        | gtf_gene_filter.run (
            runIf: {id, state -> !state.transcript_fasta}, 
            fromState: ["fasta": "fasta", "gtf": "gtf"], 
            toState: ["filtered_gtf": "filtered_gtf"],
            args: [filtered_gtf: "genome_additional_genes.gtf"] 
        )

        // prepare reference for RSEM
        | rsem_prepare_reference.run ( 
            runIf: {id, state -> !state.transcript_fasta}, 
            fromState: [
                "fasta": "fasta", 
                "gtf": "filtered_gtf" 
            ], 
            toState: ["transcript_fasta": "transcript_fasta"], 
            key: "rsem_ref",
            args: [transcript_fasta: "transcriptome.fasta"]
        )

        // chromosome size and fai index
        | getchromsize.run (
            fromState: ["fasta": "fasta"], 
            toState: [
                "fai": "fai", 
                "sizes": "sizes" ], 
            key: "chromsizes", 
            args: [ 
                fai: "genome_additional.fasta.fai", 
                sizes: "genome_additional.fasta.sizes"
            ] 
        )
        
        // untar bbsplit index, if available
        | untar.run (
            runIf: {id, state -> state.bbsplit_index}, 
            fromState: ["input": "bbsplit_index"], 
            toState: ["bbsplit_index": "output"], 
            key: "bbsplit_uncompressed",
            args: [output: "bbsplit_index"] 
        )
        
        // create bbsplit index, if not already availble
        | bbmap_bbsplit.run (
            runIf: {id, state -> !state.bbsplit_index}, 
            fromState: [ 
                "primary_ref": "fasta", 
                "bbsplit_fasta_list": "bbsplit_fasta_list" 
            ], 
            toState: ["bbsplit_index": "bbsplit_index"], 
            args: [
                "only_build_index": true, 
                bbsplit_index: "bbsplit_index"
            ], 
            key: "bbsplit_index_uncompressed" 
        )

        // Uncompress STAR index or generate from scratch if required
        | untar.run (
            runIf: {id, state -> state.star_index}, 
            fromState: ["input": "star_index"], 
            toState: ["star_index": "output"], 
            key: "star_index_uncompressed",
            args: [output: "star_index"]
        )
        
        | star_genomegenerate.run (
            runIf: {id, state -> !state.star_index}, 
            fromState: [ 
                "fasta": "fasta", 
                "gtf": "gtf" 
            ], 
            toState: ["star_index": "star_index"], 
            key: "star_index_uncompressed",
            args: [star_index: "star_index"]
        )

        // TODO: Uncompress RSEM index or generate from scratch if required

        // TODO: Uncompress HISAT2 index or generate from scratch if required

        // Uncompress Salmon index or generate from scratch if required
        | untar.run (
            runIf: {id, state -> state.salmon_index}, 
            fromState: ["input": "salmon_index"], 
            toState: ["salmon_index": "output"], 
            key: "salmon_index_uncompressed",
            args: [output: "salmon_index"]
        )

        | salmon_index.run (
            runIf: {id, state -> !state.salmon_index}, 
            fromState: [ 
                "genome_fasta": "fasta", 
                "transcriptome_fasta": "transcript_fasta" 
            ], 
            toState: ["salmon_index": "salmon_index"], 
            key: "salmon_index_uncompressed",
            args: [salmon_index: "salmon_index"] 
        )

        | setState ( 
            [ "fasta_uncompressed": "fasta", 
            "gtf_uncompressed": "gtf", 
            "transcript_fasta_uncompressed": "transcript_fasta", 
            "gene_bed_uncompressed": "gene_bed",
            "star_index_uncompressed": "star_index", 
            "salmon_index_uncompressed": "salmon_index", 
            "bbsplit_index_uncompressed": "bbsplit_index", 
            "chrom_sizes": "sizes", 
            "fai": "fai" ] )      

    emit: 
        output_ch
}
