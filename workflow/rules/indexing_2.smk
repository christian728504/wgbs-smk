rule indexing:
    output:
        gem_index="results/indexes/reference.gem",
        fasta="results/indexes/reference.fasta",
        fasta_fai="results/indexes/reference.fasta.fai",
        chromsizes="results/indexes/chromsizes.tsv",
        sam_dict="results/indexes/sam_dict.tsv",
        indexes=directory("results/indexes"),
        signal="results/indexes/.continue"
    params:
        indexes_tar=config['indexes_tar'],
        tmpdir=config['tmpdir']
    resources:
        threads=8
        cores=4
        mem_mb=16000
    container: "docker://clarity001/wgbs-smk:latest"
    log: "results/logfiles/indexes.log"
    shell:
        """
        exec >> {log} 2>&1

        OUTPUT_DIR={params.tmpdir}/results
        mkdir -p "$OUTPUT_DIR"

        curl -L {params.indexes_tar} | lz4 -dc | tar -xvf - -C $OUTPUT_DIR/.
        cp -t {output.indexes} "$OUTPUT_DIR"/indexes/*

        touch {output.signal}
        """

        # rm -rf "$OUTPUT_DIR/indexes"
