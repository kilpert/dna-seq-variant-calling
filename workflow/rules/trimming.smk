rule cutadapt:
    input:
        get_cutadapt_input,
    output:
        fastq1="results/trimmed/{sample}/{unit}_R1.fastq.gz",
        fastq2="results/trimmed/{sample}/{unit}_R2.fastq.gz",
        qc="results/trimmed/{sample}/{unit}.paired.qc.txt",
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    benchmark:
        "benchmarks/cutadapt/{sample}-{unit}.log"
    params:
        extra="",
        adapters=get_cutadapt_adapters,
    threads: 8
    resources:
        mem_mb=2048,
    wrapper:
        "v1.12.0/bio/cutadapt/pe"


rule merge_fastqs:
    input:
        get_fastqs,
    output:
        "results/merged/{sample}.{read}.fastq.gz",
    log:
        "logs/merge-fastqs/{sample}.{read}.log"
    benchmark:
        "benchmarks/merge-fastqs/{sample}.{read}.log",
    wildcard_constraints:
        read="R1|R2",
    resources:
        mem_mb=128,
    shell:
        "cat {input} > {output} 2> {log}"


rule cutadapt_pipe:
    input:
        get_cutadapt_pipe_input,
    output:
        pipe("pipe/cutadapt/{sample}/{unit}.{fq}.fastq.gz"),
    log:
        "logs/pipe-fastqs/catadapt/{sample}-{unit}.{fq}.fastq.gz.log"
    benchmark:
        "benchmarks/pipe-fastqs/catadapt/{sample}-{unit}.{fq}.fastq.gz.log",

    threads: 0  # this does not need CPU
    resources:
        mem_mb=128,
    shell:
        "cat {input} > {output} 2> {log}"