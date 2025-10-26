import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Calculate contig mapping depth and generate plot.")
parser.add_argument("--inspector", type=Path, required=True, help="")
parser.add_argument("--output", type=Path, required=True, help="")
args = parser.parse_args()

contig_length_df = pd.read_csv(args.inspector / "contig_length_info", sep="\t", header=None, names=["contig", "length"])

contig_list = []
mapped_read_length_list = []
depth_list = []

for row in contig_length_df.itertuples():
    contig = row.contig
    length = row.length
    large_file = args.inspector / "map_depth" / f"maplength_large_{contig}"
    small_file = args.inspector / "map_depth" / f"maplength_small_{contig}"
    if large_file.exists():
        maplength_file = large_file
    elif small_file.exists():
        maplength_file = small_file
    else:
        raise FileNotFoundError(f"Neither {large_file} nor {small_file} exists.")
    
    with open(maplength_file, "r") as f:
        mapped_read_length = int(f.read().strip())
        depth = mapped_read_length / length
        contig_list.append(contig)
        mapped_read_length_list.append(mapped_read_length)
        depth_list.append(depth)

args.output.mkdir(parents=True, exist_ok=True)

depth_df = pd.DataFrame({"contig": contig_list, "mapped_read_length": mapped_read_length_list, "depth": depth_list})
df = contig_length_df.merge(depth_df, on="contig")
df.to_csv(args.output / "contig_info.tsv", sep="\t", index=False)

average_depth = sum(df["mapped_read_length"]) / sum(df["length"])
with open(args.output / "average_depth.txt", "w") as f:
    f.write(f"{average_depth}\n")

fig, ax = plt.subplots()
ax.scatter(df["length"], df["depth"], alpha=0.5, color='black')
ax.axhline(average_depth, color='red', linestyle='dashed', linewidth=1, label=f'Average depth: {average_depth:.2f}')
ax.legend()
ax.set_xscale("log")
ax.set_title("Contig length vs. mapping depth")
ax.set_xlabel("Contig length")
ax.set_ylabel("Mapping depth")
plt.savefig(args.output / "contig_depth.png", dpi=300)
plt.savefig(args.output / "contig_depth.pdf", dpi=300)
plt.close()
