#!/usr/bin/env python3
from __future__ import annotations

import sys
import re
from typing import List, Optional

HEADER_RE = re.compile(r"^\s*(SW|score|\*)", re.IGNORECASE)

def parse_out_line_to_outxm(line: str) -> Optional[List[str]]:
    line = line.rstrip("\n")
    if not line.strip():
        return None
    if HEADER_RE.match(line):
        return None

    f = line.split()
    if not f:
        return None

    # RepeatMasker .out data line:
    # score div del ins query q_begin q_end (q_left) strand repeat class/family rep_begin rep_end (rep_left) ID [*]
    # -> 15 (no star) or 16 (with star)
    if len(f) not in (15, 16):
        sys.stderr.write(f"[WARN] Unexpected field count ({len(f)}), skipping: {line}\n")
        return None

    has_star = (len(f) == 16)
    star = "*" if has_star else ""

    score, div, dele, ins = f[0:4]
    query = f[4]
    q_begin, q_end, q_left = f[5:8]
    strand = f[8] # '+' or 'C'
    repeat_name = f[9]
    class_family = f[10]
    rep_begin, rep_end, rep_left = f[11:14]
    # f[14] is ID (discard)
    # f[15] is '*' if present

    repeat_full = f"{repeat_name}#{class_family}"

    # out.xm: 13 cols (+ optional '*')
    row = [
        score, div, dele, ins,
        query, q_begin, q_end, q_left,
        strand, repeat_full,
        rep_begin, rep_end, rep_left,
    ]
    if has_star:
        row.append(star)  # 14th col
    return row

def out_to_outxm(in_path: str, out_path: str) -> None:
    n_out = 0
    n_skip = 0
    with open(in_path, "r", encoding="utf-8", errors="replace") as fin, \
         open(out_path, "w", encoding="utf-8") as fout:
        for line in fin:
            row = parse_out_line_to_outxm(line)
            if row is None:
                n_skip += 1
                continue
            fout.write("\t".join(row) + "\n")
            n_out += 1
    sys.stderr.write(f"[INFO] wrote {n_out} rows (skipped {n_skip}) to {out_path}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} repeatmasker.out repeatmasker.out.xm")
    out_to_outxm(sys.argv[1], sys.argv[2])