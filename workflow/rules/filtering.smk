rule filter_by_maf:
    input:
        "results/calls/{group}.snpsift.bcf"
    output:
        "results/calls/{group}.snpsift.maf.{maf}.bcf"
    # wildcard_constraints:
    #     dataset="\d+"
    conda:
        "../envs/vembrane.yaml"
    shell:
        """vembrane filter "INFO.get('AF', 0) <= {wildcards.maf} and QUAL > 10" {input} -o {output}"""