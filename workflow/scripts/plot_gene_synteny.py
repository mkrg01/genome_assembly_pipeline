#!/usr/bin/env python3
"""Plot the same MCScanX gene anchors in gene-order and physical coordinates."""
import argparse
import csv
import json
from pathlib import Path

from gene_synteny import read_genes, write_json, write_table
from plot_self_dotplot import (
    plt, FuncFormatter, add_contig_label_axes, build_contig_offsets,
    build_boundary_segments, add_line_collection, genome_axis_label,
)


def coordinates(genes, contigs, mode):
    counts = {c['name']: c['genes'] for c in contigs}
    offsets, rows, total = build_contig_offsets([
        (c['name'], c['length'] if mode == 'mb' else max(1, counts[c['name']])) for c in contigs])
    points = {}
    for gene in genes:
        if gene['contig'] in offsets:
            local = (gene['start'] + gene['end']) / 2 if mode == 'mb' else gene['rank'] + 0.5
            points[gene['gene_id']] = offsets[gene['contig']] + local
    return points, rows, total


def draw(genes, anchors, contigs, mode, organism_name, min_length, output):
    points, rows, total = coordinates(genes, contigs, mode)
    fig, ax = plt.subplots(figsize=(180 / 25.4, 180 / 25.4), dpi=144)
    selected = [r for r in anchors if r['gene_a'] in points and r['gene_b'] in points]
    # Count an unordered gene pair once, even if it belongs to overlapping blocks.
    plotted = {}
    for row in selected:
        key = tuple(sorted((row['gene_a'], row['gene_b'])))
        plotted.setdefault(key, set()).add(row['orientation'])
    try:
        if not rows:
            ax.text(0.5, 0.5, 'No contigs passed the length filter', ha='center', va='center', transform=ax.transAxes)
            ax.set_axis_off()
        else:
            add_line_collection(ax, build_boundary_segments(rows, total), '#CBCBCB', 0.25, 0.65, 1, 'dotted')
            ax.plot([0, total], [0, total], color='#CBCBCB', linewidth=0.5, zorder=1)
            # Orientation is the MCScanX block direction, not the strand of an individual gene.
            for orientation, color, label in [('+', '#222222', 'Forward block'), ('-', '#2674AD', 'Reverse block'), ('mixed', '#9769A8', 'Both block directions')]:
                pairs = [p for p, directions in plotted.items() if (orientation == 'mixed' and len(directions) > 1) or directions == {orientation}]
                if pairs:
                    x = [points[a] for a, b in pairs] + [points[b] for a, b in pairs]
                    y = [points[b] for a, b in pairs] + [points[a] for a, b in pairs]
                    ax.scatter(x, y, s=1.5, c=color, linewidths=0, label=label, zorder=3)
            if not plotted:
                ax.text(0.5, 0.5, 'No collinear anchors on the displayed contigs', ha='center', va='center', transform=ax.transAxes, color='#555555')
            else:
                ax.legend(loc='upper left', markerscale=3, fontsize=7)
            unit = 'Mb' if mode == 'mb' else 'gene order'
            label = genome_axis_label(organism_name, min_length).replace('genome (Mb', f'genome ({unit}')
            ax.set(xlim=(0, total), ylim=(0, total), xlabel=label, ylabel=label)
            ax.set_aspect('equal', adjustable='box')
            if mode == 'mb':
                ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x / 1e6:g}'))
                ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x / 1e6:g}'))
            ax.tick_params(length=3, labelsize=7)
            add_contig_label_axes(ax, rows)
        fig.savefig(output, dpi=600)
    finally:
        plt.close(fig)
    return dict(mode=mode, contigs=len(contigs), unique_anchor_pairs=len(plotted), plotted_points=2 * len(plotted), total_axis_length=total)


def plot(prepared, blocks, output_dir, organism_name, min_long_contig_length):
    prepared, root = Path(prepared), Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    genes = read_genes(prepared / 'genes.tsv')
    meta = json.loads((prepared / 'inputs.json').read_text())
    contigs = [c for c in meta['contigs'] if c['length'] >= min_long_contig_length]
    with (Path(blocks) / 'anchors.tsv').open() as handle:
        anchors = list(csv.DictReader(handle, delimiter='\t'))
    figures = [draw(genes, anchors, contigs, mode, organism_name, min_long_contig_length, root / f'{mode}.pdf') for mode in ('gene_order', 'mb')]
    assert figures[0]['unique_anchor_pairs'] == figures[1]['unique_anchor_pairs']
    write_table(root / 'contigs.tsv', ['name', 'length', 'genes'], contigs)
    write_json(root / 'plot.json', dict(schema_version=1, min_long_contig_length=min_long_contig_length, figures=figures,
        gene_order='Ranks of all representative genes within each contig, concatenated in decreasing contig length; a gene-free contig occupies one empty slot',
        mb='Gene midpoints in the genome; concatenated in decreasing contig length',
        symmetry='Each unordered anchor pair plotted on both sides of the diagonal',
        filtering='Contig length affects display only; analysis uses all representative proteins'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('prepared', 'blocks', 'output-dir'):
        parser.add_argument('--' + key, type=Path, required=True)
    parser.add_argument('--organism-name', required=True)
    parser.add_argument('--min-long-contig-length', type=int, default=1_000_000)
    args = vars(parser.parse_args())
    if args['min_long_contig_length'] < 0:
        parser.error('Minimum contig length must be nonnegative')
    plot(**args)


if __name__ == '__main__':
    main()
