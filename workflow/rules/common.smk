import pandas as pd

samples = (
    pd.read_csv(
        config["samples"],
        sep="\t",
        dtype={"sample_name": str, "group": str},
        comment="#",
    )
    .set_index("sample_name", drop=False)
    .sort_index()
)

units = (
    pd.read_csv(
        config["units"],
        sep="\t",
        dtype={"sample_name": str, "unit_name": str},
        comment="#",
        skip_blank_lines=True,
    )
    .set_index(["sample_name", "unit_name"], drop=False)
    .sort_index()
)

groups = samples["group"].unique()

datatype = "dna"
species = config["ref"]["species"]
build = config["ref"]["build"]
release = config["ref"]["release"]
genome_name = f"genome.{datatype}.{species}.{build}.{release}"
genome_prefix = f"results/resources/{genome_name}"
genome = f"{genome_prefix}.fasta"
genome_fai = f"{genome}.fai"
genome_dict = f"{genome_prefix}.dict"
target = config.get("target", None)

wildcard_constraints:
    group="|".join(groups),
    sample="|".join(samples["sample_name"]),

def get_group_samples(group):
    return samples.loc[samples["group"] == group]["sample_name"]

def get_group(sample):
    return samples.loc[sample]["group"]

def get_cutadapt_input(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.unit]
    return expand(
        "pipe/cutadapt/{S}/{U}.{{read}}.fastq.gz".format(
            S=unit.sample_name, U=unit.unit_name
        ),
        read=["fq1", "fq2"],
    )

def get_cutadapt_adapters(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.unit]
    try:
        adapters = unit["adapters"]
        if isinstance(adapters, str):
            return adapters
        return ""
    except KeyError:
        return ""

def get_fastqs(wc):
    return expand(
        "results/trimmed/{sample}/{unit}_{read}.fastq.gz",
        unit=units.loc[wc.sample, "unit_name"],
        sample=wc.sample,
        read=wc.read,
    )

def get_cutadapt_pipe_input(wildcards):
    return get_raw_reads(wildcards.sample, wildcards.unit, wildcards.fq)


def get_raw_reads(sample, unit, fq):
    pattern = units.loc[sample].loc[unit, fq]
    if pd.isna(pattern):
        return []

    if "*" in pattern:
        files = sorted(glob.glob(units.loc[sample].loc[unit, fq]))
        if not files:
            raise ValueError(
                "No raw fastq files found for unit pattern {} (sample {}). "
                "Please check the your sample sheet.".format(unit, sample)
            )
    else:
        files = [pattern]
    return files


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample, platform=samples.loc[wildcards.sample, "platform"]
    )