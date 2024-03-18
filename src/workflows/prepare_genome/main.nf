workflow run_wf {

    take: 
        input_ch

    main: 
        output_ch = input_ch 

        // Uncompress fasta
        | gunzip.run (
            fromState: [
                "input": "fasta", 
                "versions": "versions" 
            ], 
            toState: [ 
                "fasta": "output", 
                "versions": "updated_versions" 
            ], 
            key: "gunzip_fasta",
            args: [ output: "reference_genome.fasta" ] 
        )

        // uncompress gtf
        | gunzip.run ( 
            runIf: {id, state -> state.gtf},
            fromState: [
                "input": "gtf", 
                "versions": "versions" 
            ], 
            toState: [
                "gtf": "output", 
                "versions": "updated_versions" 
            ], 
            key: "gunzip_gtf",
            args: [output: "gene_annotation.gtf"]
        )

        // uncompress gff
        | gunzip.run ( 
            runIf: {id, state -> !state.gtf && state.gff},
            fromState: [
                "input": "gff", 
                "versions": "versions" 
            ], 
            toState: [
                "gff": "output", 
                "versions": "updated_versions" 
            ], 
            key: "gunzip_gff",
            args: [output: "gene_annotation.gff"] 
        )

        // gff to gtf
        | gffread.run (
            runIf: {id, state -> !state.gtf && state.gff}, 
            fromState: [
                "input": "annotation", 
                "versions": "versions" 
            ], 
            toState: [
                "gtf": "output", 
                "versions": "updated_versions" 
            ],
            args: [output: "gene_annotation.gtf"] 
        )

        | gtf_filter.run(
            runIf: {id, state -> state.gtf && state.filter_gtf}, 
            fromState: [
                "fasta": "fasta", 
                "gtf": "gtf", 
                "versions": "versions"
            ], 
            toState: [
                "gtf": "filtered_gtf", 
                "versions": "updated_versions"
            ],
            args: [filtered_gtf: "gene_annotation.gtf"]
        )

        // uncompress additional fasta
        | gunzip.run (
            runIf: {id, state -> state.additional_fasta}, 
            fromState: [
                "input": "additional_fasta", 
                "versions": "versions" 
            ], 
            toState: [
                "additional_fasta": "output", 
                "versions": "updated_versions" 
            ], 
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
                "biotype": "biotype", 
                "versions": "versions"  
            ], 
            toState: [
                "fasta": "fasta_output", 
                "gtf": "gtf_output", 
                "versions": "updated_versions" 
            ], 
            key: "cat_additional",
            args: [
                fasta_output: "genome_additional.fasta", 
                gtf_output: "genome_additional.gtf"]  
        ) 

        // uncompress bed file
        | gunzip.run (
            runIf: {id, state -> state.gene_bed}, 
            fromState: [
                "input": "gene_bed", 
                "versions": "versions"
            ], 
            toState: [
                "gene_bed": "output", 
                "versions": "updated_versions"
            ], 
            key: "gunzip_gene_bed",
            args: [output: "genome_additional.bed"]
        )

        // gtf to bed 
        | gtf2bed.run (
            runIf: { id, state -> !state.gene_bed}, 
            fromState: [
                "gtf": "gtf", 
                "versions": "versions"
            ], 
            toState: [
                "gene_bed": "bed_output", 
                "versions": "updated_versions"
            ], 
            args: [bed_output: "genome_additional.bed"]
        ) 

        // uncompress transcript fasta
        | gunzip.run (
            runIf: {id, state -> state.transcript_fasta}, 
            fromState: [
                "input": "transcript_fasta", 
                "versions": "versions"
            ], 
            toState: [
                "transcript_fasta": "output", 
                "versions": "updated_versions"
            ], 
            key: "transcript_fasta", 
            args: [output: "transcriptome.fasta"]
        )

        // preprocess transcripts fasta if gtf is in gencode format
        | preprocess_transcripts_fasta.run (
            runIf: {id, state -> state.transcript_fasta && state.gencode}, 
            fromState: [
                "transcript_fasta": "transcript_fasta", 
                "versions": "versions"
            ], 
            toState: [
                "transcript_fasta": "output", 
                "versions": "updated_versions"
            ], 
            key: "transcript_fixed", 
            args: [output: "transcriptome.fasta"] 
        )

        // filter gtf for genes in genome
        | gtf_gene_filter.run (
            runIf: {id, state -> !state.transcript_fasta}, 
            fromState: [
                "fasta": "fasta", 
                "gtf": "gtf", 
                "versions": "versions"
            ], 
            toState: [
                "filtered_gtf": "filtered_gtf", 
                "versions": "updated_versions"
            ],
            args: [filtered_gtf: "genome_additional_genes.gtf"] 
        )

        // prepare reference for RSEM
        | rsem_prepare_reference.run ( 
            runIf: {id, state -> !state.transcript_fasta}, 
            fromState: [
                "fasta": "fasta", 
                "gtf": "filtered_gtf", 
                "versions": "versions" 
            ], 
            toState: [
                "transcript_fasta": "transcript_fasta", 
                "versions": "updated_versions"
            ], 
            key: "rsem_ref",
            args: [transcript_fasta: "transcriptome.fasta"]
        )

        // chromosome size and fai index
        | getchromsizes.run (
            fromState: [
                "fasta": "fasta", 
                "versions": "versions"
            ], 
            toState: [
                "fai": "fai", 
                "sizes": "sizes", 
                "versions": "updated_versions" 
            ], 
            key: "chromsizes", 
            args: [ 
                fai: "genome_additional.fasta.fai", 
                sizes: "genome_additional.fasta.sizes"
            ] 
        )
        
        // untar bbsplit index, if available
        | untar.run (
            runIf: {id, state -> state.bbsplit_index}, 
            fromState: [
                "input": "bbsplit_index", 
                "versions": "versions"
            ], 
            toState: [
                "bbsplit_index": "output", 
                "versions": "updated_versions"
            ], 
            key: "bbsplit_uncompressed",
            args: [output: "BBSplit_index"] 
        )
        
        // create bbsplit index, if not already availble
        | bbmap_bbsplit.run (
            runIf: {id, state -> !state.skip_bbsplit && !state.bbsplit_index}, 
            fromState: [ 
                "primary_ref": "fasta", 
                "bbsplit_fasta_list": "bbsplit_fasta_list", 
                "versions": "versions"
            ], 
            toState: [
                "bbsplit_index": "bbsplit_index", 
                "versions": "updated_versions"
            ], 
            args: [
                "only_build_index": true, 
                bbsplit_index: "BBSplit_index"
            ], 
            key: "bbsplit_index_uncompressed" 
        )

        // Uncompress STAR index or generate from scratch if required
        | untar.run (
            runIf: {id, state -> state.star_index}, 
            fromState: [
                "input": "star_index", 
                "versions": "versions"
            ], 
            toState: [
                "star_index": "output", 
                "versions": "updated_versions"
            ], 
            key: "star_index_uncompressed",
            args: [output: "STAR_index"]
        )
        
        | star_genomegenerate.run (
            runIf: {id, state -> !state.star_index}, 
            fromState: [ 
                "fasta": "fasta", 
                "gtf": "gtf", 
                "versions": "versions" 
            ], 
            toState: [
                "star_index": "star_index", 
                "versions": "updated_versions"
            ], 
            key: "star_index_uncompressed",
            args: [star_index: "STAR_index"]
        )

        // TODO: Uncompress RSEM index or generate from scratch if required

        // TODO: Uncompress HISAT2 index or generate from scratch if required

        // Uncompress Salmon index or generate from scratch if required
        | untar.run (
            runIf: {id, state -> state.salmon_index}, 
            fromState: [
                "input": "salmon_index", 
                "versions": "versions"
            ], 
            toState: [
                "salmon_index": "output", 
                "versions": "updated_versions"
            ], 
            key: "salmon_index_uncompressed",
            args: [output: "Salmon_index"]
        )

        | salmon_index.run (
            runIf: {id, state -> !state.salmon_index}, 
            fromState: [ 
                "genome_fasta": "fasta", 
                "transcriptome_fasta": "transcript_fasta", 
                "versions": "versions" 
            ], 
            toState: [
                "salmon_index": "salmon_index", 
                "versions": "updated_versions"
            ], 
            key: "salmon_index_uncompressed",
            args: [salmon_index: "Salmon_index"] 
        )

        | setState ( 
            "fasta_uncompressed": "fasta", 
            "gtf_uncompressed": "gtf", 
            "transcript_fasta_uncompressed": "transcript_fasta", 
            "gene_bed_uncompressed": "gene_bed",
            "star_index_uncompressed": "star_index", 
            "salmon_index_uncompressed": "salmon_index", 
            "bbsplit_index_uncompressed": "bbsplit_index", 
            "chrom_sizes": "sizes", 
            "fai": "fai", 
            "updated_versions": "versions"
        )      

    emit: 
        output_ch
}
