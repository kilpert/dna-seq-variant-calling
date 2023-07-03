rule table:
    input:
        "results/calls/{group}.snpsift.bcf"
    output:
        "results/tables/{group}.tsv"
    params:
        expression=lambda wc: "CHROM, POS, REF, ALT, QUAL, CSQ['SYMBOL'], CSQ['Consequence'], CSQ['IMPACT'], CSQ['Feature'], INFO['AF'], INFO['AC'], INFO['nhomalt']," + ", ".join(f"FORMAT['GT']['{s}'], FORMAT['AD']['{s}'][0], FORMAT['AD']['{s}'][1]" for s in get_group_samples(wc.group))
    conda:
        "../envs/vembrane.yaml"
    benchmark:
        "benchmarks/table/{group}.txt"
    resources:
        mem_mb=256
    shell:
        """vembrane table --overwrite-number-format GT=2 --annotation-key CSQ "{params.expression}" {input} > {output}"""


rule table_maf:
    input:
        "results/calls/{group}.snpsift.maf.{maf}.bcf"
    output:
        "results/tables/{group}.maf.{maf}.tsv"
    params:
        expression=lambda wc: "CHROM, POS, REF, ALT, QUAL, CSQ['SYMBOL'], CSQ['Consequence'], CSQ['IMPACT'], CSQ['Feature'], INFO['AF'], INFO['AC'], INFO['nhomalt']," + ", ".join(f"FORMAT['GT']['{s}'], FORMAT['AD']['{s}'][0], FORMAT['AD']['{s}'][1]" for s in get_group_samples(wc.group))
    conda:
        "../envs/vembrane.yaml"
    benchmark:
        "benchmarks/table_maf/{group}.{maf}.txt"
    resources:
        mem_mb=256
    shell:
        """vembrane table --overwrite-number-format GT=2 --annotation-key CSQ "{params.expression}" {input} > {output}"""
