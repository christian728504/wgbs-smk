rule gembs_extract:
    input:
        indexes_singal="results/indexes/.continue",
        mapping_signal="results/mapping/{barcode}/.continue",
        calls_signal="results/calls/{barcode}/.continue",
        bcf="results/calls/{barcode}/{barcode}.bcf",
        bcf_csi="results/calls/{barcode}/{barcode}.bcf.csi",
        chromsizes=rules.indexing.output.chromsizes,
    output:
        signal="results/extract/{barcode}/.continue",
        cpgfile="results/extract/{barcode}/{barcode}_cpg.txt.gz",
        noncpgfile="results/extract/{barcode}/{barcode}_non_cpg.txt.gz",
        bed_gz=[
            "results/extract/{barcode}/{barcode}_chg.bed.gz",
            "results/extract/{barcode}/{barcode}_chh.bed.gz",
            "results/extract/{barcode}/{barcode}_cpg.bed.gz"
        ],
        bigwigs=[
            "results/extract/{barcode}/{barcode}_neg.bw",
            "results/extract/{barcode}/{barcode}_pos.bw",
        ]
    params:
        tmpdir=config['tmpdir'],
        barcode=lambda w: w.barcode,
        phred_threshold=f'--threshold {config['extract']['phred_threshold']}' if config['extract']['phred_threshold'] else '',
        mode='--mode strand-specific' if config['extract']['strand_specific'] else '',
        bw_mode='--bw-mode strand-specific' if config['extract']['bw_strand_specific'] else ''
    container: "docker://clarity001/wgbs-smk:latest"
    log: "results/logfiles/extract/{barcode}.log"
    shell:
        """
        exec >> {log} 2>&1

        OUTPUT_DIR={params.tmpdir}/results/extract/{params.barcode}
        mkdir -p "$OUTPUT_DIR"

        cp {input.bcf} {input.bcf_csi} {input.chromsizes} $OUTPUT_DIR

        mextr \
            --loglevel info \
            --compress \
            --md5 \
            --cpgfile $OUTPUT_DIR/$(basename {output.cpgfile}) \
	          --noncpgfile $OUTPUT_DIR/$(basename {output.noncpgfile}) \
            --bed-methyl $OUTPUT_DIR/{params.barcode} \
            --tabix \
            --threads {resources.threads} \
            {params.phred_threshold} \
            {params.mode} \
            {params.bw_mode} \
            $OUTPUT_DIR/$(basename {input.bcf})

        rm -rf $OUTPUT_DIR/$(basename {input.chromsizes}) $OUTPUT_DIR/$(basename {input.bcf}) $OUTPUT_DIR/$(basename {input.bcf_csi})
        cp -t results/extract/{params.barcode} "$OUTPUT_DIR"/* 
        rm -rf "$OUTPUT_DIR"

        touch {output.signal}
        """

rule extract_signal:
    input:
        signal="results/extract/{barcode}/.continue",
        chromsizes=rules.indexing.output.chromsizes,
        bed_gz=[
            "results/extract/{barcode}/{barcode}_chg.bed.gz",
            "results/extract/{barcode}/{barcode}_chh.bed.gz",
            "results/extract/{barcode}/{barcode}_cpg.bed.gz"
        ],
    output:
        signal="results/extract_signal/{barcode}/.continue",
        bigwigs=[
            "results/extract_signal/{barcode}/{barcode}_cpg_neg.bw",
            "results/extract_signal/{barcode}/{barcode}_cpg_pos.bw",
            "results/extract_signal/{barcode}/{barcode}_chg_pos.bw",
            "results/extract_signal/{barcode}/{barcode}_chg_neg.bw",
            "results/extract_signal/{barcode}/{barcode}_chh_pos.bw",
            "results/extract_signal/{barcode}/{barcode}_chh_neg.bw",
        ]
    container: "docker://clarity001/wgbs-smk:latest"
    log: "results/logfiles/extract_signal/{barcode}.log"
    shell:
        """
        python workflow/rules/scripts/extract_signal.py \
            --bed-gz {input.bed_gz[0]} \
                     {input.bed_gz[1]} \
                     {input.bed_gz[2]} \
            --chromsizes {input.chromsizes} \
            --output-dir results/extract_signal/{wildcards.barcode} \
            --barcode {wildcards.barcode} \
            --log-file {log}
        """

rule get_methylc:
    input:
        signal="results/mapping/{barcode}/.continue",
        chh_bed_gz="results/extract/{barcode}/{barcode}_chh.bed.gz",
        chg_bed_gz="results/extract/{barcode}/{barcode}_chg.bed.gz",
        cpg_bed_gz="results/extract/{barcode}/{barcode}_cpg.bed.gz",
        chromsizes=rules.indexing.output.chromsizes
    output:
        signal="results/get_methylc/{barcode}/.continue",
        methylc_bed_gz="results/get_methylc/{barcode}/{barcode}_methylc.bed.gz",
    container: "docker://clarity001/wgbs-smk:latest"
    log: "results/logfiles/get_methylc/{barcode}.log"
    shell:
        """
        exec >> {log} 2>&1

        python workflow/rules/scripts/get_methylc.py \
            --barcode {wildcards.barcode} \
            --chh {input.chh_bed_gz} \
            --chg {input.chg_bed_gz} \
            --cpg {input.cpg_bed_gz} \
            --chromsizes {input.chromsizes} \
            --threads {resources.threads} \
            --outdir results/get_methylc/{wildcards.barcode} \

        touch {output.signal}
        """
