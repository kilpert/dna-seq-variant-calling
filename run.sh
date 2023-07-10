#!/usr/bin/env bash

snakemake --profile humgen_slurm --rerun-triggers mtime -j 64 -s /projects/humgen/pipelines/dna-seq-strelka.testing/workflow/Snakefile

