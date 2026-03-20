# This pipeline is adapted from and references DongAhn’s segmental duplication pipeline.

import pandas as pd
import gzip
from pathlib import Path

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

def get_asm_manifest_df(manifest_df):
    add_haps = {"H2": "hap2"}
    rows = []
    for idx, row in manifest_df.iterrows():
        if "H1" in manifest_df.columns and pd.notna(row["H1"]) and str(row["H1"]).strip():
            rows.append({"SAMPLE": row["SAMPLE"], "HAP": "hap1", "FASTA": row["H1"], "RM":row["RM"], "TRF":row["TRF"]})
        for col in add_haps:
            if col in manifest_df.columns and pd.notna(row[col]) and str(row[col]).strip():
                rows.append({"SAMPLE": row["SAMPLE"], "HAP": add_haps[col], "FASTA": row[col], "RM":row["RM"], "TRF":row["TRF"]})
    return pd.DataFrame(rows)

raw_manifest_df = pd.read_csv("manifest.tab", sep="\t", comment="#")
manifest_df = get_asm_manifest_df(raw_manifest_df)
groups = manifest_df[["SAMPLE","HAP"]].drop_duplicates().copy()
manifest_df.set_index(["SAMPLE","HAP"], inplace=True)


include: "rules/windowmasker.smk"
include: "rules/mask.smk"
include: "rules/sedef.smk"
