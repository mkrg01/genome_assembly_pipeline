from snakemake.script import snakemake
import pandas as pd
import matplotlib.pyplot as plt

def load_coverage_df(input_file):
    cols = ["contig", "start", "end", "n_feature", "n_base", "window_size", "coverage"]
    df = pd.read_csv(input_file, sep="\t", header=None, names=cols)
    global window_size_global
    window_size_global = df["window_size"].iloc[0]
    return df[["contig", "start", "end", "n_feature"]].rename(columns={"n_feature": "count"})

def load_tidk_df():
    df = pd.read_csv(snakemake.input.tidk, sep="\t")
    df = df.rename(columns={"id": "contig"})
    window_size = df["window"][0]
    global window_size_global
    window_size_global = window_size
    df["end"] = df["window"]
    df["start"] = (df["end"] - 1) // window_size * window_size
    df["count"] = df["forward_repeat_number"] + df["reverse_repeat_number"]
    return df[["contig", "start", "end", "count"]]

def load_track_df(track_cfg):
    track_id = track_cfg["id"]
    if track_id == "gene":
        return load_coverage_df(snakemake.input.gene)
    if track_id in {"LTR", "Copia", "Gypsy", "LINE", "SINE", "DNA_transposon", "satellite"}:
        return load_coverage_df(snakemake.input[track_id])
    if track_id == "tidk":
        return load_tidk_df()
    raise ValueError(f"Unsupported track id: {track_id}")

# Initialize global window size variable
window_size_global = 10000  # Default value

# Load contig information
contig_df = pd.read_csv(snakemake.input.contig, sep="\t", header=None, names=["contig", "length"])
contigs = contig_df["contig"].tolist()
contig_lengths = dict(zip(contig_df["contig"], contig_df["length"]))

# Load track data
circos_tracks = [track for track in snakemake.config["circos_plot_tracks"]]
n_tracks = len(circos_tracks)

track_data = {}
for track_cfg in circos_tracks:
    track_id = track_cfg["id"]
    track_data[track_id] = {
        "df": load_track_df(track_cfg),
        "label": track_cfg["label"],
        "color": track_cfg["color"]
    }

# Create figure with subplots for each contig and track
n_contigs = len(contigs)
fig_height = max(10, n_tracks * 1.5)
fig_width = max(4 * n_contigs, 20)
fig, axes = plt.subplots(n_tracks, n_contigs, figsize=(fig_width, fig_height), squeeze=False)

# Calculate x-axis ticks for each contig (to be shared across tracks)
from matplotlib.ticker import MaxNLocator
contig_x_ticks = {}
for contig in contigs:
    contig_len = contig_lengths[contig]
    locator = MaxNLocator(nbins=5, integer=True)
    ticks = locator.tick_values(0, contig_len)
    # Filter to keep only ticks within the range
    ticks = [tick for tick in ticks if 0 <= tick <= contig_len]
    contig_x_ticks[contig] = ticks

# Plot each track and each contig
for row_idx, track_cfg in enumerate(circos_tracks):
    track_id = track_cfg["id"]
    
    # Find y-axis max across all contigs for this track
    all_y_for_track = track_data[track_id]["df"]["count"].to_numpy()
    if len(all_y_for_track) > 0:
        y_max_global = max(all_y_for_track) * 1.1
    else:
        y_max_global = 1
    
    for col_idx, contig in enumerate(contigs):
        ax = axes[row_idx, col_idx]
        contig_len = contig_lengths[contig]
        
        # Filter data for this contig
        contig_data = track_data[track_id]["df"]
        contig_data = contig_data[contig_data["contig"] == contig].copy()
        
        if not contig_data.empty:
            # Plot the count as filled area with local contig coordinates
            x = ((contig_data["start"] + contig_data["end"]) / 2).to_numpy()
            y = contig_data["count"].to_numpy()
            ax.fill_between(x, y, color=track_cfg["color"], alpha=0.7)
        
        # Set x-axis limits to actual contig length
        ax.set_xlim(0, contig_len)
        
        # Set x-axis ticks to use shared tick positions for this contig
        ax.set_xticks(contig_x_ticks[contig])
        
        # Use global y-max for all contigs to enable comparison
        ax.set_ylim(0, y_max_global)
        
        # Set background color
        ax.set_facecolor((0.95, 0.95, 0.95))
        ax.grid(True, alpha=0.3, axis='y')
        
        # Label the top row with contig names
        if row_idx == 0:
            # Adjust font size based on contig name length and size
            name_len = len(contig)
            font_size = max(6, min(9, 80 / max(name_len, 1)))
            ax.set_title(contig, fontsize=font_size, fontweight='bold')
        
        # Y-axis label on the left side
        if col_idx == 0:
            if window_size_global >= 1_000_000 and window_size_global % 1_000_000 == 0:
                window_mb = window_size_global // 1_000_000
                y_label = f"Count per\n{int(window_mb)}-Mb window"
            elif window_size_global >= 1000 and window_size_global % 1000 == 0:
                window_kb = window_size_global // 1000
                y_label = f"Count per\n{int(window_kb)}-kb window"
            else:
                y_label = f"Count per\n{int(window_size_global)}-bp window"
            ax.set_ylabel(y_label, rotation=90, ha='center', va='center', fontsize=8, labelpad=15)
            
            # Add track label to the left of the y-axis label
            ax.text(-0.22, 0.5, track_cfg["label"], transform=ax.transAxes, 
                   fontsize=9, fontweight='bold', color=track_cfg["color"],
                   ha='center', va='center', rotation=90)
        else:
            ax.set_ylabel("")
        
        # X-axis formatting
        if row_idx == n_tracks - 1:
            # Show x-axis labels and convert to Mb
            from matplotlib.ticker import FuncFormatter
            def format_mb(x, pos):
                mb = x / 1_000_000
                if mb == int(mb):
                    return f"{int(mb)}"
                else:
                    return f"{mb:.1f}"
            ax.xaxis.set_major_formatter(FuncFormatter(format_mb))
            ax.tick_params(axis='x', labelsize=7)
            ax.set_xlabel("Position (Mb)", fontsize=8)
        else:
            ax.set_xlabel("")
            ax.set_xticklabels([])
        
        # Y-axis formatting
        ax.tick_params(axis='y', labelsize=7)
        
        # Remove spines for cleaner look
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)

plt.tight_layout()
fig.savefig(snakemake.output[0], dpi=300, bbox_inches="tight")
plt.close()