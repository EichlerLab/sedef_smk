# This pipeline is adapted from and references DongAhn’s segmental duplication pipeline.

import pandas as pd
import gzip
from pathlib import Path
import sys
import os

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

def get_asm_manifest_df(raw_manifest_df):
    hap_map = {
        "H1": "hap1",
        "H2": "hap2",
    }

    rows = []
    for _, row in raw_manifest_df.iterrows():
        sample = row["SAMPLE"]

        for col, hap in hap_map.items():
            if col in raw_manifest_df.columns and pd.notna(row[col]) and str(row[col]).strip():
                rows.append(
                    {
                        "SAMPLE": sample,
                        "HAP": hap,
                        "FASTA": row[col],
                    }
                )

    return pd.DataFrame(rows)


raw_manifest_df = pd.read_csv("manifest.tab", sep="\t", comment="#")

# sample-level manifest
sample_manifest_df = raw_manifest_df.set_index("SAMPLE", drop=False).copy()

# hap-level manifest
asm_manifest_df = get_asm_manifest_df(raw_manifest_df)
asm_manifest_df.set_index(["SAMPLE", "HAP"], inplace=True)

groups = asm_manifest_df.reset_index()[["SAMPLE", "HAP"]].drop_duplicates().copy()

rule all:
    input:
        expand("results/{sample}/final_outputs/beds/{hap}.SDs.merged.bed",
            zip,
            sample=groups["SAMPLE"],
            hap=groups["HAP"],
        ),
        expand("results/{sample}/final_outputs/bigbeds/{hap}.SDs.bb",
            zip,
            sample=groups["SAMPLE"],
            hap=groups["HAP"],
        )

include: "rules/windowmasker.smk"
include: "rules/mask.smk"
include: "rules/sedef.smk"