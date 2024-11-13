workflow run_wf {

    take: 
        input_ch

    main: 
        output_ch = input_ch 

        // Uncompress fasta
        | gunzip.run (
            fromState: [ "input": "fasta" ], 
            toState: [ "fasta": "output" ], 
            key: "gunzip_fasta",
            args: [ output: "reference_genome.fasta" ] 
        )

        // uncompress gtf
        | gunzip.run ( 
            runIf: {id, state -> state.gtf},
            fromState: [ "input": "gtf" ], 
            toState: [ "gtf": "output" ], 
            key: "gunzip_gtf",
            args: [output: "gene_annotation.gtf"]
        )

        // uncompress gff
        | gunzip.run ( 
            runIf: {id, state -> !state.gtf && state.gff},
            fromState: [ "input": "gff" ], 
            toState: [ "gff": "output" ], 
            key: "gunzip_gff",
            args: [output: "gene_annotation.gff"] 
        )

        // gff to gtf
        | gffread.run (
            runIf: {id, state -> !state.gtf && state.gff}, 
            fromState: [ 
                "input": "gff",
                "genome": "fasta" 
            ], 
            toState: [ "gtf": "outfile" ],
            args: [
                outfile: "gene_annotation.gtf", 
                gtf_output: true
            ] 
        )

        | gtf_filter.run(
            runIf: {id, state -> state.gtf && state.filter_gtf}, 
            fromState: [
                "fasta": "fasta", 
                "gtf": "gtf"
            ], 
            toState: [ "gtf": "filtered_gtf" ],
            args: [filtered_gtf: "gene_annotation.gtf"]
        )

        // uncompress additional fasta
        | gunzip.run (
            runIf: {id, state -> state.additional_fasta}, 
            fromState: [ "input": "additional_fasta" ], 
            toState: [ "additional_fasta": "output" ], 
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
            args: [
                fasta_output: "genome_additional.fasta", 
                gtf_output: "genome_additional.gtf"
            ]  
        ) 

        // uncompress bed file
        | gunzip.run (
            runIf: {id, state -> state.gene_bed}, 
            fromState: [ "input": "gene_bed" ], 
            toState: [ "gene_bed": "output" ], 
            key: "gunzip_gene_bed",
            args: [output: "genome_additional.bed"]
        )

        // gtf to bed 
        | gtf2bed.run (
            runIf: { id, state -> !state.gene_bed}, 
            fromState: [ "gtf": "gtf" ], 
            toState: [ "gene_bed": "bed_output" ], 
            args: [bed_output: "genome_additional.bed"]
        ) 

        // uncompress transcript fasta
        | gunzip.run (
            runIf: {id, state -> state.transcript_fasta}, 
            fromState: [ "input": "transcript_fasta" ], 
            toState: [ "transcript_fasta": "output" ], 
            key: "transcript_fasta", 
            args: [output: "transcriptome.fasta"]
        )

        // preprocess transcripts fasta if gtf is in gencode format
        | preprocess_transcripts_fasta.run (
            runIf: {id, state -> state.transcript_fasta && state.gencode}, 
            fromState: [ "transcript_fasta": "transcript_fasta" ], 
            toState: [ "transcript_fasta": "output" ], 
            args: [output: "transcriptome.fasta"] 
        )

        // make transcript FASTA if not provided
        | rsem_prepare_reference.run ( 
            runIf: {id, state -> !state.transcript_fasta}, 
            fromState: [
                "reference_fasta_files": "fasta", 
                "gtf": "gtf"
            ], 
            toState: [ "make_transcript_fasta_output": "output" ], 
            key: "make_transcript_fasta",\
            args: [reference_name: "genome"]
        )
        | map { id, state -> 
            def transcript_fasta = (!state.transcript_fasta) ?
                state.make_transcript_fasta_output.listFiles().find{it.name == "genome.transcripts.fa"} : 
                state.transcript_fasta
            [ id, state + [transcript_fasta: transcript_fasta] ]
        }

        // chromosome size and fai index
        | getchromsizes.run (
            fromState: [ "fasta": "fasta" ], 
            toState: [
                "fai": "fai", 
                "sizes": "sizes"
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
            fromState: [ "input": "bbsplit_index" ], 
            toState: [ "bbsplit_index": "output" ], 
            key: "untar_bbsplit_index",
            args: [output: "BBSplit_index"] 
        )

        | map {id, state -> 
            def ref = [state.fasta] + state.bbsplit_fasta_list
            [id, state + [bbsplit_ref: ref] ]
        }

        // create bbsplit index, if not already availble
        | bbmap_bbsplit.run (
            runIf: {id, state -> !state.skip_bbsplit && !state.bbsplit_index}, 
            fromState: ["ref": "bbsplit_ref"],
            toState: [ "bbsplit_index": "index" ], 
            args: [
                only_build_index: true,
                index: "BBSplit_index"
            ], 
            key: "generate_bbsplit_index"
        )

        // Uncompress STAR index or generate from scratch if required
        | untar.run (
            runIf: {id, state -> state.star_index}, 
            fromState: [ "input": "star_index" ], 
            toState: [ "star_index": "output" ], 
            key: "untar_star_index",
            args: [output: "STAR_index"]
        )
        
        // TODO: Add to viah-hub or adapt star_align_reads to enable the generateGenome runMode 
        | star_genome_generate.run (
            runIf: {id, state -> !state.star_index && !state.skip_alignment}, 
            fromState: [ 
                "genome_fasta_files": "fasta", 
                "sjdb_gtf_file": "gtf"
            ], 
            toState: [ "star_index": "index" ], 
            key: "generate_star_index",
            args: [index: "STAR_index"]
        )

        // Uncompress RSEM index or generate from scratch if required
        | untar.run (
            runIf: {id, state -> state.rsem_index}, 
            fromState: [ "input": "rsem_index" ], 
            toState: [ "rsem_index": "output" ], 
            key: "untar_rsem_index",
            args: [output: "RSEM_index"]
        )

        | rsem_prepare_reference.run ( 
            runIf: {id, state -> !state.rsem_index && state.aligner == 'star_rsem'}, 
            fromState: [
                "reference_fasta_files": "fasta", 
                "gtf": "gtf"
            ], 
            toState: [ "rsem_index": "output" ], 
            key: "generate_rsem_index",
            args: [reference_name: "genome"]
        )
        
        // TODO: Uncompress HISAT2 index or generate from scratch if required

        // Uncompress Salmon index or generate from scratch if required
        | untar.run (
            runIf: {id, state -> state.salmon_index}, 
            fromState: [ "input": "salmon_index" ], 
            toState: [ "salmon_index": "output" ], 
            key: "untar_salmon_index",
            args: [output: "Salmon_index"]
        )

        | salmon_index.run (
            runIf: {id, state -> (state.aligner == 'star_salmon' || state.pseudo_aligner == "salmon") && !state.salmon_index}, 
            fromState: [ 
                "genome": "fasta", 
                "transcripts": "transcript_fasta", 
                "kmer_len": "pseudo_aligner_kmer_size",
                "gencode": "gencode"
            ], 
            toState: [ "salmon_index": "index" ], 
            key: "generate_salmon_index",
            args: [index: "Salmon_index"] 
        )

        // Uncompress Kallisto index or generate from scratch if required
        | untar.run (
            runIf: {id, state -> state.kallisto_index}, 
            fromState: [ "input": "kallisto_index" ], 
            toState: [ "kallisto_index": "output" ], 
            key: "untar_kallisto_index",
            args: [output: "Kallisto_index"]
        )

        | kallisto_index.run(
            runIf: {id, state -> state.pseudo_aligner == "kallisto" && !state.kallisto_index}, 
            fromState: [
                "input": "transcript_fasta",
                "kmer_size": "pseudo_aligner_kmer_size"
            ],
            toState: [ "kallisto_index": "index" ],
            key: "generate_kallisto_index",
            args: [index: "Kallisto_index"]
        )

        | map { id, state -> 
          def mod_state = state.findAll { key, value -> value instanceof java.nio.file.Path && value.exists() }
          [ id, mod_state ]
        }

        | setState ( 
            "fasta_uncompressed": "fasta", 
            "gtf_uncompressed": "gtf", 
            "transcript_fasta_uncompressed": "transcript_fasta", 
            "gene_bed_uncompressed": "gene_bed",
            "star_index_uncompressed": "star_index", 
            "salmon_index_uncompressed": "salmon_index",
            "kallisto_index_uncompressed": "kallisto_index", 
            "bbsplit_index_uncompressed": "bbsplit_index", 
            "rsem_index_uncompressed": "rsem_index",
            "chrom_sizes": "sizes", 
            "fai": "fai"
        )      

    emit: 
        output_ch
}
