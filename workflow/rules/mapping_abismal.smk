
rule mapping:
    input:
        signal="results/indexes/.continue",
        abismal_index=rules.indexing.output.abismal_index,
        chromsizes=rules.indexing.output.chromsizes,
        r1=get_r1,
        r2=get_r2
    output:
        mapping=directory("results/mapping/{barcode}"),
        bam="results/mapping/{barcode}/{barcode}.bam",
        bam_csi="results/mapping/{barcode}/{barcode}.bam.csi",
        md5="results/mapping/{barcode}/{barcode}.bam.md5",
        json="results/mapping/{barcode}/{barcode}.json",
        signal="results/mapping/{barcode}/.continue"
    params:
        dataset=lambda w: dataset_lookup.get(w.barcode),
        barcode=lambda w: w.barcode,
        tmpdir=config['tmpdir']
    container: "docker://clarity001/wgbs-smk:latest"
    log: "results/logfiles/mapping/{barcode}.log"
    shell:
        """
        exec >> {log} 2>&1

        echo "$(date): Rule started"

        OUTPUT_DIR={params.tmpdir}/{output.mapping}
        
        mkdir -p "$OUTPUT_DIR"
        cp {input.abismal_index} {input.chromsizes} {input.r1} {input.r2} $OUTPUT_DIR

        echo "$(date): Copied files to /tmp"

        abismal map \
            -v \
            -t {resources.threads} \
            -i $OUTPUT_DIR/$(basename {input.abismal_index}) \
            -B -o $OUTPUT_DIR/$(basename {output.bam}) \
            $OUTPUT_DIR/$(basename {input.r1}) $OUTPUT_DIR/$(basename {input.r2}) | \
        read_filter.py $OUTPUT_DIR/$(basename {input.chromsizes}) | \
        samtools sort -o $OUTPUT_DIR/$(basename {output.bam}) \
            -T "$OUTPUT_DIR/sort" \
            --threads {resources.threads} \
            --write-index -
        
        echo "$(date): gem-mapper completed"
        
        md5sum $OUTPUT_DIR/$(basename {output.bam}) > $OUTPUT_DIR/$(basename {output.md5})

        echo "$(date): md5sum completed"

        make_average_coverage.py \
            --bamfile $OUTPUT_DIR/$(basename {output.bam}) \
            --chromsizes $OUTPUT_DIR/$(basename {input.chromsizes}) \
            --threads {resources.threads} \
            --gem_mapper_json $OUTPUT_DIR/$(basename {output.json})

        echo "$(date): average coverage completed"

        rm -rf $OUTPUT_DIR/$(basename {input.r1}) \
                $OUTPUT_DIR/$(basename {input.r2}) \
                $OUTPUT_DIR/$(basename {input.abismal_index}) \
                $OUTPUT_DIR/$(basename {input.chromsizes})
        cp -t {output.mapping} "$OUTPUT_DIR"/* 
        rm -rf "$OUTPUT_DIR"
        
        echo "$(date): Rule finished"
        touch {output.signal}
        """