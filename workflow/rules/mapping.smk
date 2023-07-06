rule map_reads:
    threads: 64
    input:
        reads=["results/merged/{sample}.R1.fastq.gz", "results/merged/{sample}.R2.fastq.gz"],
        reference=genome,
        idx=rules.bwa_index.output
    output:
        "results/mapped/{sample}.cram"
    params:
        bwa="bwa-mem2",
        extra=lambda wc: get_read_group(wc) + " -M",
        sort="samtools",  # Can be 'none' or 'samtools'.
        sort_order="coordinate",  # Can be 'coordinate' (default) or 'queryname'.
        sort_extra="",  # Extra args for samtools.
        dedup="mark",  # Can be 'none' (default), 'mark' or 'remove'.
        dedup_extra="-M",  # Extra args for samblaster.
        exceed_thread_limit=True,  # Set threads als for samtools sort / view (total used CPU may exceed threads!)
        embed_ref=True,  # Embed reference when writing cram.
    log:
        "logs/bwa_meme/{sample}.log",
    benchmark:
        "benchmarks/map_reads/{sample}.txt"
    resources:
        mem_mb=150000
    wrapper:
        "v1.25.0/bio/bwa-memx/mem"


def gatk_baserecalibratorspar_extra():
    if target:
        return f"-L {target}"
    return ""

rule gatk_baserecalibratorspark:
    input:
        bam="results/mapped/{sample}.cram",
        bai="results/mapped/{sample}.cram.crai",
        ref=genome,
        dict=genome_dict,
        known="/projects/humgen/science/resources/gnomad.genomes.v3.1.2.sites.genome.vcf.gz",
    output:
        recal_table="results/recal_table/{sample}.grp",
    params:
        extra=gatk_baserecalibratorspar_extra(),
        java_opts="",  # optional
        #spark_runner="",  # optional, local by default
        #spark_v1.25.0="",  # optional
        #spark_extra="", # optional
    log:
        "logs/gatk/baserecalibrator/{sample}.log",
    benchmark:
        "benchmarks/gatk/baserecalibrator/{sample}.txt"
    resources:
        mem_mb=10000,
    threads: 8
    wrapper:
        "v2.0.0/bio/gatk/baserecalibratorspark"


rule gatk_applybqsr_spark:
    input:
        bam="results/mapped/{sample}.cram",
        bai="results/mapped/{sample}.cram.crai",
        ref=genome,
        dict=genome_dict,
        recal_table="results/recal_table/{sample}.grp",
    output:
        bam="results/alignment/{sample}.cram",
    params:
        extra="--create-output-bam-index false",  # optional,
        #spark_master="local[32]",
        #spark_runner="",  # optional, local by default
        #spark_v1.25.0="",  # optional
        #spark_extra="", # optional
        embed_ref=True,  # embed reference in cram output
        exceed_thread_limit=True,  # samtools is also parallized and thread limit is not guaranteed anymore
    log:
        "logs/gatk_applybqsr_spark/{sample}.log",
    benchmark:
        "benchmarks/gatk_applybqsr_spark/{sample}.txt"
    resources:
        mem_mb=90000,
    threads: 32
    wrapper:
        "v2.0.0/bio/gatk/applybqsrspark"


rule link_alignment_index:
    input:
        "{x}.crai"
    output:
        "{x}.cram.crai"
    resources:
        mem_mb=128,
    shell:
        "ln {input} {output}"


rule samtools_index:
    input:
        "{x}.cram",
    output:
        "{x}.crai",
    params:
        extra="",  # optional params string
    log:
        "logs/samtools_index/{x}.log"
    benchmark:
        "benchmarks/samtools_index/{x}.log"
    resources:
        mem_mb=1024,
    threads: 8  # This value - 1 will be sent to -@
    wrapper:
        "v1.25.0/bio/samtools/index"
