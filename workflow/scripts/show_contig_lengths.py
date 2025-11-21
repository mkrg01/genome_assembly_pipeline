import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser(description="Plot contig lengths.")
parser.add_argument("--input", type=Path, required=True, help="")
parser.add_argument("--output", type=Path, required=True, help="")
args = parser.parse_args()

contig_df = pd.read_csv(args.input, sep="\t")

# Summarise contigs with less than 1Mbp length, renaming contig names to "<1Mbp (n contigs)"
short_contigs = contig_df[contig_df["length"] < 1_000_000]
n_short_contigs = short_contigs.shape[0]
if n_short_contigs > 0:
    short_contig_summary = pd.DataFrame({
        "#name": [f"<1Mbp ({n_short_contigs} contigs)"],
        "length": [short_contigs["length"].sum()]
    })
    contig_df = pd.concat([
        contig_df[contig_df["length"] >= 1_000_000],
        short_contig_summary
    ], ignore_index=True)

fig, ax = plt.subplots()
ax.bar(contig_df["#name"], contig_df["length"], color='black')
ax.set_title("Contig length distribution")
ax.set_xlabel("Contigs")
x = np.arange(len(contig_df))
ax.set_xticks(x)
ax.set_xticklabels(contig_df["#name"], rotation=90)
ax.set_ylabel("Lengths (bp)")
plt.tight_layout()
plt.savefig(args.output, dpi=300)
plt.close()
