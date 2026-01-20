rule generate_chrome_ctgs:
    input:
        indexes_signal="results/indexes/.continue",
        chromsizes=rules.indexing.output.chromsizes
    output:
        ctg_beds=temp(expand("results/calls/{{barcode}}/{{barcode}}_{contig}_ctgs.bed", contig=CONTIGS))
    params:
        barcode=lambda w: w.barcode,
        tmpdir=config['tmpdir']
    container: "docker://clarity001/wgbs-smk:latest"
    log: "results/logfiles/calls/{barcode}_ctgs.log"
    shell:
        """
        exec >> {log} 2>&1

        OUTPUT_DIR={params.tmpdir}/results/calls/{params.barcode}
        mkdir -p "$OUTPUT_DIR"
        cp {input.chromsizes} $OUTPUT_DIR

        python workflow/rules/scripts/make_ctg_beds.py \
            --chromsizes $OUTPUT_DIR/$(basename {input.chromsizes}) \
            --output_dir $OUTPUT_DIR \
            --barcode {params.barcode}

        rm -rf $OUTPUT_DIR/$(basename {input.chromsizes})
        cp -t results/calls/{params.barcode} "$OUTPUT_DIR"/{params.barcode}_*_ctgs.bed
        rm -rf "$OUTPUT_DIR"
        
        echo "Generated contig bed files!"
        """

rule calling:
    input:
        bam="results/mapping/{barcode}/{barcode}.bam",
        bam_csi="results/mapping/{barcode}/{barcode}.bam.csi",
        ctg_bed="results/calls/{barcode}/{barcode}_{contig}_ctgs.bed",
        fasta=rules.indexing.output.fasta,
        fasta_fai=rules.indexing.output.fasta_fai,
        chromsizes=rules.indexing.output.chromsizes
    output:
        bcf=temp("results/calls/{barcode}/{barcode}_{contig}.bcf"),
        report="results/calls/{barcode}/{barcode}_{contig}.json"
    params:
        barcode=lambda w: w.barcode,
        contig=lambda w: w.contig,
        tmpdir=config['tmpdir'],
        out_type=config['calling']['output_type'],
        conversion=config['calling']['conversion'],
        ltrim=config['calling']['left_trim'],
        rtrim=config['calling']['right_trim'],
        ref_bias=config['calling']['reference_bias'],
        mapq=config['calling']['mapq_threshold'],
        qual=config['calling']['qual_threshold']
    container: "docker://clarity001/wgbs-smk:latest"
    log: "results/logfiles/calls/{barcode}_{contig}.log"
    shell:
        """
        exec >> {log} 2>&1

        OUTPUT_DIR={params.tmpdir}/results/calls/{params.barcode}/{params.contig}
        CTGS_DIR={params.tmpdir}/results/calls/{params.barcode}
        mkdir -p "$OUTPUT_DIR"
        mkdir -p "$CTGS_DIR"
        cp {input.bam} {input.bam_csi} {input.chromsizes} {input.fasta} {input.fasta_fai} $OUTPUT_DIR
        cp {input.ctg_bed} $CTGS_DIR

        bs_call \
            --loglevel info \
            --output $OUTPUT_DIR/$(basename {output.bcf}) \
            --output-type {params.out_type} \
            --reference $OUTPUT_DIR/$(basename {input.fasta}) \
            --sample {params.barcode} \
            --contig-include $OUTPUT_DIR/$(basename {input.chromsizes}) \
            --report-file $OUTPUT_DIR/$(basename {output.report}) \
            --contig-bed $CTGS_DIR/$(basename {input.ctg_bed}) \
            --threads {resources.threads} \
            --conversion {params.conversion} \
            --left-trim {params.ltrim} \
            --right-trim {params.rtrim} \
            --reference-bias {params.ref_bias} \
            --mapq-threshold {params.mapq} \
            --bq-threshold {params.qual} \
            $OUTPUT_DIR/$(basename {input.bam})

        cp $OUTPUT_DIR/$(basename {output.bcf}) $OUTPUT_DIR/$(basename {output.report}) results/calls/{params.barcode}
        rm -rf $OUTPUT_DIR/$(basename {input.fasta}) $OUTPUT_DIR/$(basename {input.fasta}) $OUTPUT_DIR/$(basename {input.bam}) $OUTPUT_DIR/$(basename {input.bam_csi})
        rm -rf $CTGS_DIR/$(basename {input.ctg_bed})

        echo "Generated BCF file for {params.contig}"
        """

rule merge_chroms:
    input:
        chrom_bcfs=get_chrom_bcfs
    output:
        signal="results/calls/{barcode}/.continue",
        bcf="results/calls/{barcode}/{barcode}.bcf",
        bcf_csi="results/calls/{barcode}/{barcode}.bcf.csi",
        bcf_md5="results/calls/{barcode}/{barcode}.bcf.md5"
    params:
        barcode=lambda w: w.barcode,
        tmpdir=config['tmpdir'],
        out_type=config['calling']['output_type']
    container: "docker://clarity001/wgbs-smk:latest"
    log: "results/logfiles/calls/{barcode}_merge.log"
    shell:
        """
        exec >> {log} 2>&1

        OUTPUT_DIR={params.tmpdir}/results/calls/{params.barcode}
        mkdir -p "$OUTPUT_DIR"
        cp {input.chrom_bcfs} $OUTPUT_DIR

        bcftools concat \
            --output $OUTPUT_DIR/$(basename {output.bcf}) \
            --output-type {params.out_type} \
            --naive \
            --threads {resources.threads} \
           $OUTPUT_DIR/{params.barcode}_*.bcf

        bcftools index --threads {resources.threads} $OUTPUT_DIR/$(basename {output.bcf})

        md5sum $OUTPUT_DIR/$(basename {output.bcf}) > $OUTPUT_DIR/$(basename {output.bcf_md5})

        cp -t results/calls/{params.barcode} $OUTPUT_DIR/$(basename {output.bcf}) $OUTPUT_DIR/$(basename {output.bcf_csi}) $OUTPUT_DIR/$(basename {output.bcf_md5})
        rm -rf "$OUTPUT_DIR"

        echo "Done!"
        touch {output.signal}
        """
