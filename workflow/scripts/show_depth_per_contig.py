import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Extract target contig read depth and generate plot.")
parser.add_argument("--mapping_tsv", type=Path, required=True, help="")
parser.add_argument("--contig_names", type=Path, required=True, help="")
parser.add_argument("--outdir", type=Path, required=True, help="")
parser.add_argument("--prefix", type=Path, required=True, help="")
args = parser.parse_args()

with open(args.contig_names, "r") as f:
    contig_name_list = [line.strip() for line in f.readlines()]
contig_df = pd.read_csv(args.mapping_tsv, sep="\t")
contig_df = contig_df[contig_df["#rname"].isin(contig_name_list)]
contig_df["length"] = contig_df["endpos"] - contig_df["startpos"] + 1

args.outdir.mkdir(parents=True, exist_ok=True)

contig_df.to_csv(args.outdir / f"{args.prefix}_contig_info.tsv", sep="\t", index=False)

rough_average_depth = sum(contig_df["meandepth"] * contig_df["length"]) / sum(contig_df["length"])
with open(args.outdir / f"{args.prefix}_average_depth.txt", "w") as f:
    f.write(f"{rough_average_depth}\n")

fig, ax = plt.subplots()
ax.scatter(contig_df["length"], contig_df["meandepth"], alpha=0.5, color='black')
ax.axhline(rough_average_depth, color='red', linestyle='dashed', linewidth=1, label=f'Average depth: {rough_average_depth:.0f}')
ax.legend()
ax.set_xscale("log")
ax.set_title("Contig length vs. read depth")
ax.set_xlabel("Contig length")
ax.set_ylabel("Read depth")
plt.savefig(args.outdir / f"{args.prefix}_contig_depth.pdf", dpi=300)
plt.close()
