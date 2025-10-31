import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Plot GC content per contig.")
parser.add_argument("--input", type=Path, required=True, help="")
parser.add_argument("--output", type=Path, required=True, help="")
args = parser.parse_args()

contig_df = pd.read_csv(args.input, sep="\t")

fig, ax = plt.subplots()
ax.scatter(contig_df["length"], contig_df["GC"], alpha=0.5, color='black')
ax.set_xscale("log")
ax.set_title("Contig length vs. GC content")
ax.set_xlabel("Contig length")
ax.set_ylabel("GC content")
plt.savefig(args.output, dpi=300)
plt.close()
