rule annotate_snpsift:
    input:
        call="results/calls/{group}.annotated.bcf",
        database="/projects/humgen/science/resources/gnomad.genomes.v3.1.2.sites.genome.vcf.gz"
    output:
        call="results/calls/{group}.snpsift.bcf"
    log:
        "logs/snpsift/{group}.annotate.log"
    benchmark:
        "benchmarks/snpsift/{group}.annotate.txt"
    resources:
        mem_mb=8096
    threads: 3
    wrapper:
        "v1.31.1/bio/snpsift/annotate"
        # "v1.25.0/bio/snpsift/annotate"


#vembrane filter "not CSQ['Consequence']=='intergenic_variant' and (not INFO['AF'] or INFO['AF']<0.001) and QUAL>0" variants.snpsift.bcf --annotation-key CSQ | vembrane table --annotation-key "CSQ" "CHROM, POS, ALT, REF
rule annotate_vep:
    input:
        calls="results/calls/{group}.norm.bcf",
        cache="results/resources/vep/cache",
        plugins="results/resources/vep/plugins",
    output:
        calls="results/calls/{group}.annotated.bcf",
        stats="results/calls/{group}.vep.html",
    params:
        plugins="",
        extra=""
    log:
        "logs/vep/{group}.annotate.log"
    benchmark:
        "benchmarks/vep/{group}.annotate.txt"
    resources:
        mem_mb=8096
    threads: 8
    wrapper:
        "v1.24.0/bio/vep/annotate"


rule fix_header:
    input:
        variants="results/calls/{group}.vcf.gz",
    output:
        new_header=temp("results/new_header/{group}.txt")
    conda:
        "../envs/normalize.yaml"
    benchmark:
        "benchmarks/fix_header/{group}.txt"
    resources:
        mem_mb=256
    shell:
        """bcftools view -h {input.variants} | sed 's/##FORMAT=<ID=AD,Number=./##FORMAT=<ID=AD,Number=R/g' | sed 's/##FORMAT=<ID=ADF,Number=./##FORMAT=<ID=ADF,Number=R/g' | sed 's/##FORMAT=<ID=ADR,Number=./##FORMAT=<ID=ADR,Number=R/g' > {output}"""


rule normalize:
    threads:
        8
    input:
        variants="results/calls/{group}.vcf.gz",
        reference=genome,
        new_header="results/new_header/{group}.txt"
    output:
        "results/calls/{group}.norm.bcf",
    log:
        "logs/normalize/{group}.log"
    benchmark:
        "benchmarks/normalize/{group}.txt"
    conda:
        "../envs/normalize.yaml"
    resources:
        mem_mb=2048
    shell:
        "(bcftools reheader --threads {threads} {input.variants} -h {input.new_header} | bcftools norm --do-not-normalize --multiallelics -any --threads {threads} | vcfallelicprimitives -k -g | vcfstreamsort | vt normalize -n -r {input.reference} - | awk -F'\\t' '!_[$1,$2,$4,$5]++' | bcftools view --threads {threads} -Ob > {output}) 2> {log}"


rule normalize_sv:
    threads:
        8
    input:
        variants="results/calls_sv/{group}.bcf",
        reference=genome,
    output:
        "results/calls_sv/{group}.norm.bcf"
    log:
        "logs/normalize_sv/{group}.log"
    benchmark:
        "benchmarks/normalize_sv/{group}.txt"
    resources:
        mem_mb=2048
    conda:
        "../envs/normalize.yaml"
    shell:
        "(bcftools norm --do-not-normalize --multiallelics -any --threads {threads} {input.variants} | vcfallelicprimitives -k -g | vcfstreamsort | vt normalize -n -r {input.reference} - | awk -F'\\t' '!_[$1,$2,$4,$5]++' | bcftools view --threads {threads} -Ob > {output}) 2> {log}"
