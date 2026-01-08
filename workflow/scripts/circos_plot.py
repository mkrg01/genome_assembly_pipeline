from snakemake.script import snakemake
import pandas as pd
from pycirclize import Circos
from matplotlib.lines import Line2D
from collections import OrderedDict

def load_gene_coverage_df():
    cols = ["contig", "start", "end", "n_feature", "n_base", "window_size", "coverage"]
    df = pd.read_csv(snakemake.input.gene, sep="\t", header=None, names=cols)
    return df[["contig", "start", "end", "coverage"]]

def load_repeat_coverage_df(repeat_class):
    cols = ["contig", "start", "end", "n_feature", "n_base", "window_size", "coverage"]
    df = pd.read_csv(snakemake.input[repeat_class], sep="\t", header=None, names=cols)
    return df[["contig", "start", "end", "coverage"]]

def load_tidk_df():
    df = pd.read_csv(snakemake.input.tidk, sep="\t")
    df = df.rename(columns={"id": "contig"})
    window_size = df["window"][0]
    df["end"] = df["window"]
    df["start"] = (df["end"] - 1) // window_size * window_size
    df["coverage"] = (df["forward_repeat_number"] + df["reverse_repeat_number"]) / window_size
    return df[["contig", "start", "end", "coverage"]]

def load_track_df(track_cfg):
    track_id = track_cfg["id"]
    if track_id == "gene":
        return load_gene_coverage_df()
    if track_id in {"LTR", "Copia", "Gypsy", "LINE", "SINE", "DNA_transposon", "satellite"}:
        return load_repeat_coverage_df(track_id)
    if track_id == "tidk":
        return load_tidk_df()
    raise ValueError(f"Unsupported track id: {track_id}")

def track_height(idx, n_tracks, outer_r=95, inner_r=10):
    height = (outer_r - inner_r) / n_tracks
    r2 = outer_r - idx * height
    r1 = r2 - height * 0.9
    return r1, r2

def add_hist_track_for_df(circos, coverage_df, track_cfg, idx, n_tracks, gap_name="__gap__"):
    r1, r2 = track_height(idx, n_tracks)
    for sector in circos.sectors:
        track = sector.add_track((r1, r2))
        if sector.name == gap_name:
            track.axis(fc="none", ec="none")
            track.text(str(idx + 1), x=sector.size * 0.5, r=(r1 + r2) / 2, size=8, orientation="horizontal")
            continue
        track.axis(fc=(0.95, 0.95, 0.95, 0.6), ec="none")
        track.grid()
        coverage_contig_df = coverage_df[coverage_df["contig"] == sector.name]
        x = ((coverage_contig_df["start"] + coverage_contig_df["end"]) / 2).to_numpy()
        y = coverage_contig_df["coverage"].to_numpy()
        track.fill_between(x, y, color=track_cfg["color"])
        if idx == 0:
            track.xticks_by_interval(interval=5_000_000, outer=True, show_bottom_line=False, show_endlabel=False, label_formatter=lambda x: "", tick_length=0.5)
            track.xticks_by_interval(interval=10_000_000, outer=True, show_bottom_line=False, show_endlabel=False, label_formatter=lambda x: f"{int(x/1_000_000)} Mb", label_size=5, tick_length=1, label_orientation="vertical")

contig_df = pd.read_csv(snakemake.input.contig, sep="\t", header=None, names=["contig", "length"])

gap_name = "__gap__"
gap_frac = 0.01
total_len = int(contig_df["length"].sum())
gap_len = max(1, int(total_len * gap_frac))

base_sectors = OrderedDict(zip(contig_df["contig"], contig_df["length"]))
sectors = OrderedDict(base_sectors)
sectors[gap_name] = gap_len

circos = Circos(sectors, space=1)
for sector in circos.sectors:
    if sector.name == gap_name:
        continue
    sector.text(sector.name, size=8, orientation="vertical")

circos_tracks = [track for track in snakemake.config["circos_plot"]]
n_tracks = len(circos_tracks)
for idx, track_cfg in enumerate(circos_tracks):
    coverage_df = load_track_df(track_cfg)
    add_hist_track_for_df(circos, coverage_df, track_cfg, idx, n_tracks, gap_name=gap_name)

fig = circos.plotfig()
legend_handles = [Line2D([], [], linestyle="none") for _ in circos_tracks]
legend_labels = [f"{idx + 1}. {track_cfg['label']}" for idx, track_cfg in enumerate(circos_tracks)]
leg = fig.legend(handles=legend_handles, labels=legend_labels, loc="lower left", fontsize=8, frameon=False, handlelength=0)
for text, track_cfg in zip(leg.get_texts(), circos_tracks):
    text.set_color(track_cfg["color"])
    text.set_fontweight("bold")

fig.savefig(snakemake.output[0], dpi=300, bbox_inches="tight")