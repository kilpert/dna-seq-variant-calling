rule call:
    input:
        bam=lambda wc: expand("results/alignment/{sample}.cram", sample=get_group_samples(wc.group)),
        bam_index=lambda wc: expand("results/alignment/{sample}.cram.crai", sample=get_group_samples(wc.group)),
        fasta=genome,
        fasta_index=genome_fai,
    output:
        variants="results/calls/{group}.vcf.gz",
        variants_index="results/calls/{group}.vcf.gz.tbi",
    params:
        sample_genomes=lambda wc: expand("results/gvcf/{sample}.vcf.gz", sample=get_group_samples(wc.group)),
        sample_genomes_indices=lambda wc: expand("results/gvcf/{sample}.vcf.gz.tbi", sample=get_group_samples(wc.group)),
        rm_existing=True,
    log:
        "logs/strelka/{group}.log"
    benchmark:
        "benchmarks/strelka/{group}.txt"
    resources:
        mem_mb=20000
    threads:
        61
    wrapper:
        "file:///projects/humgen/user-projects/cschroeder/snakemake-wrappers/bio/strelka/germline/"
        # "v1.24.0/bio/strelka/germline"


rule touch_gvcf:
    input:
        lambda wc: f"results/calls/{get_group(wc.sample)}.vcf.gz"
    output:
        "results/gvcf/{sample}.vcf.gz"
    resources:
        mem_mb=128
    shell:
        "touch {output}"


rule structural_call:
    input:
        samples=lambda wc: expand("results/alignment/{sample}.cram", sample=get_group_samples(wc.group)),
        index=lambda wc: expand("results/alignment/{sample}.cram.crai", sample=get_group_samples(wc.group)),
        ref=genome,
        ref_index=genome_fai,
    output:
        vcf="results/calls_sv/{group}.bcf",
        idx="results/calls_sv/{group}.bcf.csi",
        cand_indel_vcf="results/sv_calls/{group}.cand_indel.vcf.gz",
        cand_indel_idx="results/sv_calls/{group}.cand_indel.vcf.gz.tbi",
        cand_sv_vcf="results/sv_calls/{group}.cand_sv.vcf.gz",
        cand_sv_idx="results/sv_calls/{group}.cand_sv.vcf.gz.tbi",
    # params:
    #     extra_cfg="--config config/manta.ini",  # optional
    log:
        "logs/manta/{group}.log",
    benchmark:
        "benchmarks/manta/{group}.txt"
    resources:
        mem_mb=20000
    threads:
        61
    wrapper:
        "v1.24.0/bio/manta"
