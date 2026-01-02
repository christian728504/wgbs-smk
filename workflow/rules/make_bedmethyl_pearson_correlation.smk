rule make_bedmethyl_pearson_correlation:
    input:
        bedmethyls=make_bedmethyl_pearson_correlation
    output:
        csv="results/{experiment}_pearson_correlation_qc.csv"
    container: "docker://clarity001/wgbs-smk:latest"
    log: "results/logfiles/rule_make_bedmethyl_pearson_correlation/{experiment}.log"
    shell:
        """
        exec >> {log} 2>&1

        python3 workflow/rules/scripts/make_bedmethyl_pearson_correlation_v2.py --bedmethyls {input.bedmethyls} --output-csv {output.csv}
        """
        
