#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

INPUT = Path("data/examples/example-differential-proteins.csv")
OUTDIR = Path("outputs/tables")
OUTDIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT)

required = {"protein_id", "gene_symbol", "log2FC", "pvalue", "padj"}
missing = required - set(df.columns)

if missing:
    raise ValueError(f"Missing required columns: {missing}")

sig = df[df["padj"] < 0.05].copy()
up = sig[sig["log2FC"] > 0].sort_values(["padj", "log2FC"], ascending=[True, False])
down = sig[sig["log2FC"] < 0].sort_values(["padj", "log2FC"], ascending=[True, True])

df.to_csv(OUTDIR / "all-differential-proteins.csv", index=False)
sig.to_csv(OUTDIR / "significant-dep.csv", index=False)
up.to_csv(OUTDIR / "top-upregulated-proteins.csv", index=False)
down.to_csv(OUTDIR / "top-downregulated-proteins.csv", index=False)

sig["gene_symbol"].dropna().drop_duplicates().to_csv(
    OUTDIR / "enrichment-input-list.txt",
    index=False,
    header=False
)

sig["protein_id"].dropna().drop_duplicates().to_csv(
    OUTDIR / "string-input-list.txt",
    index=False,
    header=False
)

print("Differential protein filtering complete.")
print(f"All proteins: {len(df)}")
print(f"Significant DEP: {len(sig)}")
print(f"Upregulated: {len(up)}")
print(f"Downregulated: {len(down)}")
