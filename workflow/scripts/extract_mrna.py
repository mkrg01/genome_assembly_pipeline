import argparse
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Extract genomic sequences corresponding to mRNA features in a GFF3 file. "
            "The output preserves introns because each mRNA interval is sliced directly "
            "from the genome FASTA."
        )
    )
    parser.add_argument("--gff3", type=Path, required=True, help="Input GFF3 file.")
    parser.add_argument("--fasta", type=Path, required=True, help="Genome FASTA file.")
    parser.add_argument("--output", type=Path, required=True, help="Output FASTA file.")
    return parser.parse_args()


def parse_attributes(attribute_text):
    attributes = {}
    for field in attribute_text.strip().split(";"):
        if not field or "=" not in field:
            continue
        key, value = field.split("=", 1)
        attributes[key] = value
    return attributes


def main():
    args = parse_args()
    genome = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))

    records = []
    with args.gff3.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "mRNA":
                continue

            seqid, _, _, start, end, score, strand, phase, attributes_text = fields
            attributes = parse_attributes(attributes_text)
            record_id = attributes.get("ID", f"{seqid}:{start}-{end}")

            if seqid not in genome:
                raise KeyError(f"{seqid} was not found in {args.fasta}")

            start_idx = int(start) - 1
            end_idx = int(end)
            sequence = genome[seqid].seq[start_idx:end_idx]

            if strand == "-":
                sequence = sequence.reverse_complement()
            elif strand != "+":
                raise ValueError(f"Unsupported strand {strand!r} for {record_id}")

            description = f"{seqid}:{start}-{end} strand={strand}"
            records.append(SeqRecord(sequence, id=record_id, description=""))

    with args.output.open("w") as handle:
        SeqIO.write(records, handle, "fasta")


if __name__ == "__main__":
    main()
