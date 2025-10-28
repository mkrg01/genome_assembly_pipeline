import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

parser = argparse.ArgumentParser(description="Detect contigs from organelles.")
parser.add_argument("--input", type=Path, required=True, help="")
parser.add_argument("--outdir", type=Path, required=True, help="")
parser.add_argument("--prefix", type=str, required=True, help="")
parser.add_argument("--min_identity", type=float, default=0.95, help="")
parser.add_argument("--min_coverage", type=float, default=0.90, help="")
args = parser.parse_args()

args.outdir.mkdir(exist_ok=True, parents=True)

mapping_df = pd.read_csv(args.input, sep="\t", header=None, usecols=range(12), names=["query_name", "query_length", "query_start", "query_end", "strand", "target_name", "target_length", "target_start", "target_end", "num_matches", "alignment_block_length", "mapping_quality"])

# Plot identity per alignment
mapping_df["identity"] = mapping_df["num_matches"] / mapping_df["alignment_block_length"]
fig, ax = plt.subplots()
ax.hist(mapping_df["identity"].clip(upper=0.999), bins=np.arange(0, 1.01, 0.01), color='grey')
ax.axvline(args.min_identity, color='red', linestyle='dashed', linewidth=1, label=f'Threshold: {args.min_identity}')
ax.legend()
ax.set_xlabel("Identity")
ax.set_ylabel("Number of alignments")
ax.set_title("Identities per alignment")
plt.savefig(args.outdir / f"{args.prefix}_mapped_alignment_identity.pdf", dpi=300)

# Calculate coverage per contig
position_dict = defaultdict(lambda: defaultdict(list))
query_length_dict = {}
for row in mapping_df.itertuples():
    if row.identity < args.min_identity or pd.isna(row.identity):
        continue
    position_dict[row.query_name][row.target_name].append((row.query_start, row.query_end))
    query_length_dict[row.query_name] = row.query_length

coverage_list = []
for query_name, pos_dict in position_dict.items():
    target_with_max_coverage = None
    max_covered_length = -1
    for target_name, positions in pos_dict.items():
        covered_length = 0
        positions.sort()
        current_start, current_end = positions[0]
        for start, end in positions[1:]:
            if start <= current_end:
                current_end = max(current_end, end)
            else:
                covered_length += current_end - current_start
                current_start, current_end = start, end
        covered_length += current_end - current_start
        if covered_length > max_covered_length:
            max_covered_length = covered_length
            target_with_max_coverage = target_name
    coverage_list.append((query_name, query_length_dict[query_name], target_with_max_coverage, max_covered_length, max_covered_length / query_length_dict[query_name]))

coverage_df = pd.DataFrame(coverage_list, columns=["query_name", "query_length", "best_target", "covered_length", "coverage"])

# Plot coverage per contig
fig, ax = plt.subplots()
ax.hist(coverage_df["coverage"].clip(upper=0.999), bins=np.arange(0, 1.01, 0.01), color='grey')
ax.axvline(args.min_coverage, color='red', linestyle='dashed', linewidth=1, label=f'Threshold: {args.min_coverage}')
ax.legend()
ax.set_xlabel("Coverage")
ax.set_ylabel("Number of contigs")
ax.set_title("Coverage per mapped contig")
plt.savefig(args.outdir / f"{args.prefix}_mapped_contig_coverage.pdf", dpi=300)

# Select organelle contigs
contigs_from_organelle_df = coverage_df[coverage_df["coverage"] >= args.min_coverage]
contigs_from_organelle_df.to_csv(args.outdir / f"{args.prefix}_organelle_contig_info.tsv", sep="\t", index=False)

# Plot organelle contig lengths
fig, ax = plt.subplots()
ax.hist(contigs_from_organelle_df["query_length"], bins=50, color='grey')
ax.set_xlabel("Contig length")
ax.set_ylabel("Number of contigs")
ax.set_title("Lengths of organelle-derived contigs")
plt.savefig(args.outdir / f"{args.prefix}_organelle_contig_length.pdf", dpi=300)

# Plot organelle contig targets
fig, ax = plt.subplots()
best_target_counts = contigs_from_organelle_df["best_target"].value_counts()
bars = ax.bar(best_target_counts.index, best_target_counts.values, color='grey')
ax.bar_label(bars, labels=best_target_counts.values, padding=1)
ax.set_xlabel("Mapped organelle")
ax.set_ylabel("Number of contigs")
ax.set_title("Organelle-derived contigs")
plt.savefig(args.outdir / f"{args.prefix}_organelle_contig_target.pdf", dpi=300)

# Save organelle contig names
with open(args.outdir / f"{args.prefix}_organelle_contigs.txt", "w") as out_f:
    for query_name in contigs_from_organelle_df["query_name"]:
        out_f.write(f"{query_name}\n")
