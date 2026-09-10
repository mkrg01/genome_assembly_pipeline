"""Select one longest CDS per GFF3 parent gene, independent of transcript naming."""

import argparse
from pathlib import Path

from filter_gff_by_fasta_ids import iter_gff_records, open_text, parse_attributes, split_parents


def select_longest(cds, gff3, output):
    parents = {}
    for _, _, fields in iter_gff_records(gff3):
        if fields is None or fields[2] != "mRNA":
            continue
        attrs = parse_attributes(fields[8])
        tx = attrs.get("ID")
        genes = split_parents(attrs.get("Parent", ""))
        if not tx or tx in parents or len(genes) != 1:
            raise ValueError("GFF3 must have unique mRNA IDs with one parent gene each")
        parents[tx] = genes[0]

    longest, seen = {}, set()

    def consider(header, sequence):
        tx = header.split()[0]
        if tx in seen or tx not in parents or not sequence:
            raise ValueError(f"Duplicate, unmapped or empty CDS record: {tx}")
        seen.add(tx)
        gene = parents[tx]
        # Ties retain the first record in the input FASTA.
        if gene not in longest or len(sequence) > len(longest[gene][1]):
            longest[gene] = (header, sequence)

    header, sequence = None, []
    with open_text(cds) as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    consider(header, "".join(sequence))
                header, sequence = line[1:].strip(), []
                if not header:
                    raise ValueError("Empty FASTA header")
            elif line.strip():
                if header is None:
                    raise ValueError("Sequence before FASTA header")
                sequence.append(line.strip())
    if header is not None:
        consider(header, "".join(sequence))
    if not seen or seen != set(parents):
        raise ValueError("CDS FASTA IDs must match the GFF3 mRNA IDs")
    with Path(output).open("w") as handle:
        for header, sequence in longest.values():
            handle.write(f">{header}\n{sequence}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cds", type=Path, required=True)
    parser.add_argument("--gff3", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    select_longest(args.cds, args.gff3, args.output)


if __name__ == "__main__":
    main()
