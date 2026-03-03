rule split:
    input:
        r1=get_r1,
        r2=get_r2,
    output:
        dir=directory("results/split/{barcode}"),
        fastqs=expand(
            "results/split/{{barcode}}/MOHD_{{barcode}}_{read}.part_{part}.fastq.gz",
            part=SPLIT,
            read=["R1", "R2"],
        ),
        signal="results/split/{barcode}/.continue",
    container: "docker://clarity001/wgbs-smk:latest"
    log: os.path.join(workflow.basedir, "results/logfiles/mapping/{barcode}.log")
    shell:
        """
        exec >> {log} 2>&1

        echo "$(date): Rule started"
  
        seqkit split2 -1 {input.r1} -2 {input.r2} -O {output.dir} -p 2 -j {resources.threads} -f
        
        echo "$(date): Rule finished"
        touch {output.signal}
        """

rule mapping:
    input:
        gem_index=rules.gem_indexer.output.gem_index,
        r1="results/split/{barcode}/MOHD_{barcode}_R1.part_{part}.fastq.gz",
        r2="results/split/{barcode}/MOHD_{barcode}_R2.part_{part}.fastq.gz",
    output:
        bam="results/mapping/{barcode}/{barcode}.{part}.bam",
        json="results/mapping/{barcode}/{barcode}.{part}.json",
        signal="results/mapping/{barcode}/.{part}.continue",
    params:
        dataset=lambda w: dataset_lookup.get(w.barcode).strip(),
        barcode=lambda w: w.barcode.strip(),
    container: "docker://clarity001/wgbs-smk:latest"
    log: os.path.join(workflow.basedir, "results/logfiles/mapping/{barcode}.{part}.log")
    shell:
        """
        exec >> {log} 2>&1
        echo "$(date): Rule started"

        gem-mapper \
            -t {resources.threads} \
            -I {input.gem_index} \
            -1 {input.r1} \
            -2 {input.r2} \
            -p \
            --report-file={output.json} \
            -r "@RG\tID:{params.dataset}\tSM:\tBC:{params.barcode}\tPU:{params.dataset}" > {output.bam}
        
        echo "$(date): Rule finished"
        touch {output.signal}
        """

rule merge:
    input:
        bams=expand(
            "results/mapping/{{barcode}}/{{barcode}}.{part}.bam",
            part=SPLIT,
        )
    output:
        bam="results/mapping/{barcode}/{barcode}.bam",
    container: "docker://clarity001/wgbs-smk:latest"
    log: os.path.join(workflow.basedir, "results/logfiles/merge/{barcode}.log")
    shell:
        """
        exec >> {log} 2>&1
        echo "$(date): Rule started"

        samtools merge \
            -@ {resources.threads} \
            {output.bam} \
            {input.bams} 

        echo "$(date): Rule finished"
        """

rule postprocess:
    input:
        bam=rules.merge.output.bam,
        jsons=expand(
            "results/mapping/{{barcode}}/{{barcode}}.{part}.json",
            part=SPLIT,
        ),
        chromsizes=rules.indexing.output.chromsizes,
    output:
        bam="results/postprocess/{barcode}/{barcode}.bam",
        bam_csi="results/postprocess/{barcode}/{barcode}.bam.csi",
        md5="results/postprocess/{barcode}/{barcode}.bam.md5",
        signal="results/postprocess/{barcode}/.continue",
    container: "docker://clarity001/wgbs-smk:latest"
    log: "results/logfiles/postprocess/{barcode}.log"
    shell:
        """
        exec >> {log} 2>&1
        echo "$(date): Rule started"
        
        python workflow/rules/scripts/read_filter.py -i {input.bam} {input.chromsizes} | \
        samtools sort -o {output.bam} \
            -T /tmp/{wildcards.barcode} \
            --threads {resources.threads} \
            --write-index -
        echo "$(date): gem-mapper completed"
        
        md5sum {output.bam} > {output.md5}
        echo "$(date): md5sum completed"

        python workflow/rules/scripts/make_average_coverage.py \
            --bamfile {output.bam} \
            --chromsizes {input.chromsizes} \
            --threads {resources.threads} \
            --gem_mapper_jsons {input.jsons}
        echo "$(date): average coverage completed"

        echo "$(date): Rule finished"
        touch {output.signal}
        """

rule get_coverage:
    input:
        bam=rules.postprocess.output.bam,
        chromsizes=rules.indexing.output.chromsizes,
        signal=rules.postprocess.output.signal,
    output:
        signal="results/get_coverage/{barcode}/.continue",
        coverage_bedgraphs=[ "results/get_coverage/{barcode}/{barcode}_coverage_neg.bg", "results/get_coverage/{barcode}/{barcode}_coverage_pos.bg", ],
        coverage_bigwigs=[
            "results/get_coverage/{barcode}/{barcode}_coverage_neg.bw",
            "results/get_coverage/{barcode}/{barcode}_coverage_pos.bw",
        ]
    container: "docker://clarity001/wgbs-smk:latest"
    log: "results/logfiles/get_coverage/{barcode}.log"
    shell:
        """
        exec >> {log} 2>&1
        echo "$(date): Rule started"

        echo "$(date --iso=seconds): Running bedtools genomecov"
        parallel -j 2 bedtools genomecov -ibam {input.bam} -bg -strand {{1}} '>' {{2}} ::: - + :::+ {output.coverage_bedgraphs}

        echo "$(date --iso=seconds): Running bedGraphToBigWig"
        parallel -j 2 bedGraphToBigWig {{1}} {input.chromsizes} {{2}} ::: {output.coverage_bedgraphs} :::+ {output.coverage_bigwigs}

        echo "$(date): Rule finished"
        touch {output.signal}
        """
