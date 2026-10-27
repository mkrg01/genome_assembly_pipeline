#!/usr/bin/env python3
"""Prepare representative genes, run DIAMOND/MCScanX, and normalize self-synteny.

Public coordinates are zero-based, half-open. MCScanX inputs use safe internal
IDs and one-based inclusive coordinates; original identifiers remain in tables.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
from urllib.parse import unquote

GENE_FIELDS = ['tool_id', 'tool_contig', 'gene_id', 'protein_id', 'contig', 'start', 'end', 'strand', 'rank']
ANCHOR_FIELDS = ['block_id', 'gene_a', 'gene_b', 'contig_a', 'contig_b', 'start_a', 'end_a', 'start_b', 'end_b', 'orientation', 'evalue']
BLOCK_FIELDS = ['block_id', 'contig_a', 'contig_b', 'orientation', 'anchors', 'start_a', 'end_a', 'start_b', 'end_b', 'score', 'evalue']

# Shared by workflow CLI entrypoints and programmatic use.
PROTEIN_EVALUE = 1e-10
MIN_PROTEIN_COVERAGE = 0.3
MAX_TARGETS = 20

# Metadata for the pinned MCScanX 1.0.0 package; no parameter overrides are passed.
MCSCANX_DEFAULTS = dict(min_anchors=5, max_gaps=25, block_evalue=1e-5, overlap_window=5)


def open_text(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else Path(path).open()


def fasta_records(path):
    seen, name, sequence = set(), None, []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(sequence)
                fields = line[1:].split()
                if not fields or fields[0] in seen:
                    raise ValueError(f'Missing or duplicate FASTA ID in {path}: {line.strip()}')
                name, sequence = fields[0], []
                seen.add(name)
            elif line.strip():
                if name is None:
                    raise ValueError(f'Sequence before FASTA header in {path}')
                sequence.append(line.strip())
        if name is not None:
            yield name, ''.join(sequence)


def genome_lengths(path):
    # Streaming: never retain a chromosome sequence in memory.
    lengths, current = {}, None
    with open_text(path) as handle:
        for line in handle:
            if line.startswith('>'):
                fields = line[1:].split()
                if not fields or fields[0] in lengths:
                    raise ValueError(f'Missing or duplicate genome FASTA ID in {path}')
                current = fields[0]
                lengths[current] = 0
            elif line.strip():
                if current is None:
                    raise ValueError(f'Sequence before FASTA header in {path}')
                lengths[current] += len(line.strip())
    if not lengths or min(lengths.values()) < 1:
        raise ValueError('Genome FASTA must contain nonempty, uniquely named sequences')
    return lengths


def write_table(path, fields, rows):
    with Path(path).open('w') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def read_genes(path):
    with Path(path).open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    for row in rows:
        for field in ('start', 'end', 'rank'):
            row[field] = int(row[field])
    return rows


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, sort_keys=True) + '\n')


def prepare(genome, proteins, gff3, output_dir):
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    lengths = genome_lengths(genome)
    features = {}
    with open_text(gff3) as handle:
        for line_no, line in enumerate(handle, 1):
            if line.startswith('##FASTA'):
                break
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) != 9:
                raise ValueError(f'Expected nine GFF3 columns at line {line_no}')
            if f[2] not in ('gene', 'mRNA', 'transcript'):
                continue
            attrs = dict(item.split('=', 1) for item in f[8].split(';') if '=' in item)
            name = unquote(attrs.get('ID', ''))
            if not name or name in features:
                raise ValueError(f'Missing or duplicate GFF3 feature ID at line {line_no}')
            start, end = int(f[3]) - 1, int(f[4])
            contig = unquote(f[0])
            if contig not in lengths or not 0 <= start < end <= lengths[contig] or f[6] not in ('+', '-'):
                raise ValueError(f'Invalid contig, coordinates, or strand at GFF3 line {line_no}')
            features[name] = dict(kind=f[2], contig=contig, start=start, end=end, strand=f[6],
                                  parents=[unquote(p) for p in attrs.get('Parent', '').split(',') if p])
    rows, seqs, represented = [], {}, set()
    for protein_id, sequence in fasta_records(proteins):
        if protein_id not in features:
            raise ValueError(f'Protein ID {protein_id} has no matching gene/transcript in GFF3')
        feature = features[protein_id]
        if feature['kind'] == 'gene':
            gene_id = protein_id
        else:
            parents = feature['parents']
            if len(parents) != 1 or parents[0] not in features or features[parents[0]]['kind'] != 'gene':
                raise ValueError(f'Transcript {protein_id} must have exactly one existing gene parent')
            gene_id = parents[0]
            gene = features[gene_id]
            if (gene['contig'], gene['strand']) != (feature['contig'], feature['strand']) or not gene['start'] <= feature['start'] < feature['end'] <= gene['end']:
                raise ValueError(f'Transcript {protein_id} is outside its gene parent')
        if gene_id in represented:
            raise ValueError(f'Multiple proteins for gene {gene_id}; supply one representative per gene')
        represented.add(gene_id)
        sequence = sequence.upper().rstrip('*')
        if not sequence or not re.fullmatch('[A-Z]+', sequence):
            raise ValueError(f'Invalid or empty protein sequence: {protein_id}')
        gene = features[gene_id]
        rows.append(dict(gene_id=gene_id, protein_id=protein_id, **{k: gene[k] for k in ('contig', 'start', 'end', 'strand')}))
        seqs[gene_id] = sequence
    if not rows:
        raise ValueError('No representative proteins were supplied')
    contigs = sorted(lengths, key=lambda name: (-lengths[name], name))
    order = {name: i for i, name in enumerate(contigs)}
    rows.sort(key=lambda r: (order[r['contig']], r['start'], r['end'], r['gene_id']))
    ranks = {}
    with (root / 'proteins.fa').open('w') as aa, (root / 'data.gff').open('w') as gff:
        for i, row in enumerate(rows, 1):
            row['tool_id'] = f'sg{i:09d}'
            row['tool_contig'] = f"sc{order[row['contig']]:09d}"
            row['rank'] = ranks.get(row['contig'], 0)
            ranks[row['contig']] = row['rank'] + 1
            aa.write(f">{row['tool_id']}\n{seqs[row['gene_id']]}\n")
            gff.write(f"sc{order[row['contig']]:09d}\t{row['tool_id']}\t{row['start'] + 1}\t{row['end']}\n")
    write_table(root / 'genes.tsv', GENE_FIELDS, rows)
    write_json(root / 'inputs.json', dict(schema_version=1, coordinate_system='0-based half-open',
        inputs={key: dict(path=str(Path(path).resolve()), bytes=Path(path).stat().st_size, mtime_ns=Path(path).stat().st_mtime_ns) for key, path in [('genome', genome), ('proteins', proteins), ('gff3', gff3)]},
        genes=len(rows), annotation_genes=sum(f['kind'] == 'gene' for f in features.values()),
        contigs=[dict(name=c, length=lengths[c], genes=ranks.get(c, 0)) for c in contigs]))


def run(command, cwd=None):
    subprocess.run(list(map(str, command)), cwd=cwd, check=True)


def search(proteins, output, threads, evalue=PROTEIN_EVALUE):
    output = Path(output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    database = output.parent / 'proteins'
    version = subprocess.check_output(['diamond', 'version'], text=True).strip()
    run(['diamond', 'makedb', '--in', proteins, '--db', database, '--threads', threads])
    command = ['diamond', 'blastp', '--query', str(Path(proteins).resolve()), '--db', str(database),
        '--out', str(output), '--threads', str(threads), '--ultra-sensitive', '--evalue', str(evalue),
        '--max-target-seqs', '0', '--max-hsps', '1', '--outfmt', '6',
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen']
    run(command)
    Path(str(database) + '.dmnd').unlink()
    write_json(output.parent / 'search.json', dict(tool=version, command=command, self_hits_removed_at='block filtering before top-N'))


def filter_hits(raw, genes, output, evalue=PROTEIN_EVALUE, min_coverage=MIN_PROTEIN_COVERAGE, max_targets=MAX_TARGETS):
    """Best HSP per directed pair, then per-query top N nonself hits; reciprocal union."""
    ids = {row['tool_id'] for row in genes}
    best, counts = {}, dict(raw_records=0, self_hits=0, failed_filters=0)
    with Path(raw).open() as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) != 14 or f[0] not in ids or f[1] not in ids:
                raise ValueError(f'Invalid DIAMOND row at line {line_no}')
            counts['raw_records'] += 1
            ev, score = float(f[10]), float(f[11])
            qlen, slen = int(f[12]), int(f[13])
            coords = [int(f[i]) for i in (6, 7, 8, 9)]
            if not math.isfinite(ev) or ev < 0 or not math.isfinite(score) or score < 0 or min(qlen, slen) < 1 or not 1 <= coords[0] <= coords[1] <= qlen or not 1 <= coords[2] <= coords[3] <= slen:
                raise ValueError(f'Invalid DIAMOND values at line {line_no}')
            if f[0] == f[1]:
                counts['self_hits'] += 1
                continue
            if ev > evalue or min((coords[1]-coords[0]+1)/qlen, (coords[3]-coords[2]+1)/slen) < min_coverage:
                counts['failed_filters'] += 1
                continue
            bucket = best.setdefault(f[0], {})
            previous = bucket.get(f[1])
            if previous is None or (-score, ev) < (-float(previous[11]), float(previous[10])):
                bucket[f[1]] = f
                if len(bucket) > max_targets:
                    worst = max(bucket, key=lambda sid: (-float(bucket[sid][11]), float(bucket[sid][10]), sid))
                    del bucket[worst]
    pairs = {}
    for query in sorted(best):
        hits = sorted(best[query].values(), key=lambda f: (-float(f[11]), float(f[10]), f[1]))
        for f in hits:
            key = tuple(sorted(f[:2]))
            if key not in pairs or (-float(f[11]), float(f[10])) < (-float(pairs[key][11]), float(pairs[key][10])):
                pairs[key] = f
    with Path(output).open('w') as handle:
        for key in sorted(pairs):
            handle.write('\t'.join(pairs[key][:12]) + '\n')
    return dict(counts, homology_pairs=len(pairs), evalue=evalue, min_protein_coverage=min_coverage,
                max_targets=max_targets, hit_selection='per-query top N nonself hits, reciprocal union; no reciprocal-best-hit restriction')


def normalize_collinearity(path, genes):
    lookup = {r['tool_id']: r for r in genes}
    header = re.compile(r'^## Alignment (\d+): score=([\d.eE+\-]+) e_value=([\d.eE+\-]+) N=(\d+) (\S+)&(\S+) (plus|minus)$')
    blocks, anchors, current = [], [], None
    seen = set()
    with Path(path).open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('## Alignment'):
                m = header.fullmatch(line)
                if not m:
                    raise ValueError(f'Invalid MCScanX block header: {line}')
                number, score, ev, count, ca, cb, orient = m.groups()
                if number in seen:
                    raise ValueError(f'Duplicate MCScanX block {number}')
                seen.add(number)
                current = dict(block_id=number, score=float(score), evalue=float(ev), orientation='+' if orient == 'plus' else '-', expected=int(count), members=[], tool_contigs=(ca, cb))
                blocks.append(current)
            elif line.startswith('#'):
                continue
            else:
                f = line.split()
                # MCScanX writes "0-  0: gene1 gene2 evalue" (4 or 5 fields).
                if current is None or len(f) < 4 or f[-3] not in lookup or f[-2] not in lookup:
                    raise ValueError(f'Invalid MCScanX anchor: {line}')
                a, b = lookup[f[-3]], lookup[f[-2]]
                if 'tool_contig' in a and (a['tool_contig'], b['tool_contig']) != current['tool_contigs']:
                    raise ValueError('MCScanX anchor contigs disagree with the block header')
                if a['gene_id'] == b['gene_id']:
                    raise ValueError('MCScanX output contains a self-gene anchor')
                value = float(f[-1])
                if not math.isfinite(value) or value < 0:
                    raise ValueError('Invalid MCScanX anchor evalue')
                row = dict(block_id=current['block_id'], gene_a=a['gene_id'], gene_b=b['gene_id'], contig_a=a['contig'], contig_b=b['contig'], start_a=a['start'], end_a=a['end'], start_b=b['start'], end_b=b['end'], orientation=current['orientation'], evalue=value)
                current['members'].append(row)
                anchors.append(row)
    normalized = []
    for b in blocks:
        members = b['members']
        if len(members) != b['expected'] or not members or len({(r['contig_a'], r['contig_b']) for r in members}) != 1:
            raise ValueError(f"Inconsistent MCScanX block {b['block_id']}")
        normalized.append(dict(block_id=b['block_id'], contig_a=members[0]['contig_a'], contig_b=members[0]['contig_b'], orientation=b['orientation'], anchors=len(members), start_a=min(r['start_a'] for r in members), end_a=max(r['end_a'] for r in members), start_b=min(r['start_b'] for r in members), end_b=max(r['end_b'] for r in members), score=b['score'], evalue=b['evalue']))
    return normalized, anchors


def blocks(prepared, raw, output_dir, evalue=PROTEIN_EVALUE, min_coverage=MIN_PROTEIN_COVERAGE, max_targets=MAX_TARGETS):
    prepared, root = Path(prepared).resolve(), Path(output_dir).resolve()
    root.mkdir(parents=True, exist_ok=True)
    genes = read_genes(prepared / 'genes.tsv')
    shutil.copyfile(prepared / 'data.gff', root / 'data.gff')
    stats = filter_hits(raw, genes, root / 'data.blast', evalue, min_coverage, max_targets)
    command = ['MCScanX', 'data', '-a']
    binary = Path(shutil.which('MCScanX')).resolve()
    fingerprint = hashlib.sha256(binary.read_bytes()).hexdigest()
    if stats['homology_pairs']:
        run(command, cwd=root)
    else:
        (root / 'data.collinearity').write_text('# No nonself protein pairs passed filters.\n')
    block_rows, anchors = normalize_collinearity(root / 'data.collinearity', genes)
    write_table(root / 'blocks.tsv', BLOCK_FIELDS, block_rows)
    write_table(root / 'anchors.tsv', ANCHOR_FIELDS, anchors)
    lookup = {g['gene_id']: g for g in genes}
    with (root / 'anchors.bedpe').open('w') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        for row in anchors:
            a, b = lookup[row['gene_a']], lookup[row['gene_b']]
            writer.writerow([a['contig'], a['start'], a['end'], b['contig'], b['start'], b['end'],
                             f"block_{row['block_id']}:{a['gene_id']}:{b['gene_id']}", 0, a['strand'], b['strand']])
    anchor_genes = {r[key] for r in anchors for key in ('gene_a', 'gene_b')}
    stats.update(anchor_genes=len(anchor_genes), anchor_gene_fraction=len(anchor_genes) / len(genes),
                 intercontig_blocks=sum(b['contig_a'] != b['contig_b'] for b in block_rows),
                 intracontig_blocks=sum(b['contig_a'] == b['contig_b'] for b in block_rows),
                 forward_blocks=sum(b['orientation'] == '+' for b in block_rows),
                 reverse_blocks=sum(b['orientation'] == '-' for b in block_rows))
    stats.update(MCSCANX_DEFAULTS)
    stats.update(schema_version=1, coordinate_system='0-based half-open', blocks=len(block_rows), anchors=len(anchors), unique_anchor_pairs=len({tuple(sorted((r['gene_a'], r['gene_b']))) for r in anchors}), genes=len(genes),
                 diamond=json.loads((Path(raw).parent / 'search.json').read_text()), mcscanx=dict(environment_package_version='1.0.0', source_tag='v1.0.0', executable=str(binary), sha256=fingerprint, command=command), interpretation='Gene collinearity; not nucleotide alignment coverage or a ploidy estimate')
    write_json(root / 'summary.json', stats)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    prep = sub.add_parser('prepare')
    for arg in ('genome', 'proteins', 'gff3', 'output-dir'):
        prep.add_argument('--' + arg, type=Path, required=True)
    srch = sub.add_parser('search')
    srch.add_argument('--proteins', type=Path, required=True)
    srch.add_argument('--output', type=Path, required=True)
    srch.add_argument('--threads', type=int, default=1)
    srch.add_argument('--evalue', type=float, default=PROTEIN_EVALUE)
    blk = sub.add_parser('blocks')
    for arg in ('prepared', 'raw', 'output-dir'):
        blk.add_argument('--' + arg, type=Path, required=True)
    blk.add_argument('--evalue', type=float, default=PROTEIN_EVALUE)
    blk.add_argument('--min-coverage', type=float, default=MIN_PROTEIN_COVERAGE)
    blk.add_argument('--max-targets', type=int, default=MAX_TARGETS)
    args = vars(parser.parse_args())
    action = args.pop('action')
    if 'evalue' in args and (not math.isfinite(args['evalue']) or args['evalue'] <= 0):
        parser.error('--evalue must be positive and finite')
    if action == 'blocks' and (not 0 <= args['min_coverage'] <= 1 or args['max_targets'] < 1):
        parser.error('Invalid coverage or target count')
    if action == 'search' and args['threads'] < 1:
        parser.error('--threads must be positive')
    globals()[action](**args)


if __name__ == '__main__':
    main()
