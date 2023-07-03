rule filter_by_maf:
    input:
        "results/calls/{group}.snpsift.bcf"
    output:
        "results/calls/{group}.snpsift.maf.{maf}.bcf"
    # wildcard_constraints:
    #     dataset="\d+"
    conda:
        "../envs/vembrane.yaml"
    benchmark:
        "benchmarks/filter_by_maf/{group}.{maf}.txt"
    resources:
        mem_mb=256
    shell:
        """vembrane filter "INFO.get('AF', 0) <= {wildcards.maf} and QUAL > 10" {input} -o {output}"""