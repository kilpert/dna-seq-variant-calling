rule map_reads:
    input:
        reads=["results/merged/{sample}.R1.fastq.gz", "results/merged/{sample}.R2.fastq.gz"],
        reference=genome,
        idx=rules.bwa_index.output
    output:
        "results/mapped/{sample}.cram"
    log:
        "logs/bwa_meme/{sample}.log",
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
    threads: 40
    wrapper:
        "v1.25.0/bio/bwa-memx/mem"
        #"file:///projects/humgen/science/depienne/dna-seq-strelka/workflow/wrapper/bwa-meme"


rule gatk_baserecalibratorspark:
    input:
        bam="results/mapped/{sample}.cram",
        bai="results/mapped/{sample}.cram.crai",
        ref=genome,
        dict=genome_dict,
        known="/projects/humgen/science/resources/gnomad.genomes.v3.1.2.sites.genome.vcf.gz",
    output:
        recal_table="results/recal_table/{sample}.grp",
    log:
        "logs/gatk/baserecalibrator/{sample}.log",
    params:
        extra=f"-L /projects/humgen/science/resources/twist2.bed",  # optional
        java_opts="",  # optional
        #spark_runner="",  # optional, local by default
        #spark_v1.25.0="",  # optional
        #spark_extra="", # optional
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
    log:
        "logs/gatk_applybqsr_spark/{sample}.log",
    params:
        extra="--create-output-bam-index false",  # optional,
        #spark_master="local[32]",
        #spark_runner="",  # optional, local by default
        #spark_v1.25.0="",  # optional
        #spark_extra="", # optional
        embed_ref=True,  # embed reference in cram output
        exceed_thread_limit=True,  # samtools is also parallized and thread limit is not guaranteed anymore
    resources:
        mem_mb=50000,
    threads: 32
    wrapper:
        "v2.0.0/bio/gatk/applybqsrspark"


rule link_alignment_index:
    input:
        "{x}.crai"
    output:
        "{x}.cram.crai"
    shell:
        "ln {input} {output}"


rule samtools_index:
    input:
        "{x}.cram",
    output:
        "{x}.crai",
    log:
        "logs/samtools_index/{x}.log",
    params:
        extra="",  # optional params string
    threads: 8  # This value - 1 will be sent to -@
    wrapper:
        "v1.25.0/bio/samtools/index"
