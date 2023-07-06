
# def samtools_stats_input(wc):
#     ret = {"bam": f"alignment/{wc.sample}.cram"}
#     

samtools_stats_input = {"bam": "results/alignment/{sample}.cram"}
if target:
    samtools_stats_input["bed"]=target

rule cram_to_bam:
    input:
        "results/alignment/{sample}.cram",
    output:
        bam=pipe("results/alignment/{sample}.bam"),
    log:
        "logs/cram_to_bam/{sample}.log"
    benchmark:
        "benchmarks/cram_to_bam/{sample}.log"
    conda:
        "../envs/samtools.yaml"
    resources:
        mem_mb=1024
    threads: 0
    shell:
        "samtools view -@ 64 {input} -O bam > {output}"


rule verify_bam_id:
    threads: 64
    input:
        bam="results/alignment/{sample}.cram",
        ref=genome,
    output:
        selfsm="results/qc/verify_bam_id/{sample}.selfSM",
        ancestry="results/qc/verify_bam_id/{sample}.ancestry",
    params:
        genome_build="38",
    log:
        "logs/verifybamid2/a{sample}.log",
    benchmark:
        "benchmarks/verifybamid2/{sample}.txt",
    wrapper:
        "v2.0.0/bio/verifybamid/verifybamid2"


rule qualimap:
    input:
        bam="results/alignment/{sample}.bam",
    output:
        outdir=directory("results/qc/qualimap/{sample}"),
    params:
        target=target
    log:
        "logs/qualimap/bamqc/{sample}.log"
    benchmark:
        "benchmarks/qualimap/bamqc/{sample}.txt"
    conda:
        "../envs/qualimap.yaml"
    threads:
        4
    resources:
        mem_gb=10,
    wrapper:
        "file:///projects/humgen/user-projects/cschroeder/snakemake-wrappers/bio/qualimap/bamqc"


rule samtools_idxstats:
    input:
        bam="results/alignment/{sample}.cram",
        idx="results/alignment/{sample}.cram.crai",
    output:
        "results/qc/samtools_idxstats/{sample}.idxstats",
    log:
        "logs/samtools_idxstats/{sample}.log"
    benchmark:
        "benchmarks/samtools_idxstats/{sample}.log"
    params:
        extra="" # optional params string
    resources:
        mem_mb=1024,
    wrapper:
        "v2.0.0/bio/samtools/idxstats"


rule samtools_stats:
    input:
        **samtools_stats_input,
    output:
        "results/qc/samtools_stats/{sample}.txt",
    params:
        extra="",  # Optional: extra arguments.
    log:
        "logs/samtools_stats/{sample}.log",
    benchmark:
        "benchmarks/samtools_stats/{sample}.log",
    resources:
        mem_mb=1024,
    wrapper:
        "v2.0.0/bio/samtools/stats"


rule multiqc_dir:
    input:
        expand("results/qc/qualimap/{sample}", sample=samples.sample_name),
        expand("results/qc/samtools_stats/{sample}.txt", sample=samples.sample_name),
        expand("results/qc/samtools_idxstats/{sample}.idxstats", sample=samples.sample_name),
        expand("results/calls/{group}.vep.html", group=groups),
        # expand("logs/bwa_meme/{sample}.log", sample=samples.sample_name),
    output:
        "results/qc/multiqc.html"
    params:
        extra="" # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log"
    benchmark:
        "benchmarks/multiqc.log"
    resources:
        mem_mb=10000,
    wrapper:
        "v2.0.0/bio/multiqc"