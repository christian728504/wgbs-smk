rule indexing:
    input:
        fasta=config['reference']
    output:
        fasta="results/indexes/reference.fasta",
        fasta_fai="results/indexes/reference.fasta.fai",
        chromsizes="results/indexes/chromsizes.tsv",
        sam_dict="results/indexes/sam_dict.tsv",
        indexes=directory("results/indexes"),
        signal="results/indexes/.continue"
    params:
        assembly=lambda w: f"-a {config.get('assembly')}" if config.get('assembly') else '',
        species=lambda w: f"-s {config.get('species')}" if config.get('species') else '',
        tmpdir=config['tmpdir']
    container: "docker://clarity001/wgbs-smk:latest"
    log: "results/logfiles/indexes.log"
    shell:
        """
        exec >> {log} 2>&1
        OUTPUT_DIR={params.tmpdir}/results/indexes
        
        mkdir -p "$OUTPUT_DIR"
        cp {input.fasta} "$OUTPUT_DIR"
        FASTA=$OUTPUT_DIR/$(basename {input.fasta})
        if [[ "$FASTA" == *.gz ]]; then
            unpigz -c "$FASTA" > $OUTPUT_DIR/$(basename {output.fasta})
        else
            mv "$FASTA" $OUTPUT_DIR/$(basename {output.fasta})
        fi
        rm -rf "$FASTA"

        samtools faidx -@ {resources.threads} $OUTPUT_DIR/$(basename {output.fasta})
        cut -f1,2 $OUTPUT_DIR/$(basename {output.fasta}).fai > $OUTPUT_DIR/$(basename {output.chromsizes})
        samtools dict \
            -o $OUTPUT_DIR/$(basename {output.sam_dict}) \
            {params.assembly} \
            {params.species} \
            $OUTPUT_DIR/$(basename {output.fasta})

        cp -t {output.indexes} "$OUTPUT_DIR"/* 
        rm -rf "$OUTPUT_DIR"

        echo "Done!"
        touch {output.signal}
        """

rule gem_indexer:
    input:
        fasta=rules.indexing.output.fasta,
    output:
        gem_index="results/indexes/reference.gem",
        signal="results/indexes/gem_indexer.continue",
    params:
        sampling_rate=config['index']['sampling_rate'],
    shadow: "copy-minimal"
    container: "docker://clarity001/gem3-mapper:latest"
    log: os.path.join(workflow.basedir, "results/logfiles/gem_indexer.log") # must be abspath
    shell:
        """
        exec >> {log} 2>&1

        gem-indexer \
            -i {input.fasta} \
            -o $(dirname {output.gem_index})/reference \
            --bisulfite-index \
            --text-sampling-rate {params.sampling_rate} \
            --threads {resources.threads}

        echo "Done!"
        touch {output.signal}
        """


