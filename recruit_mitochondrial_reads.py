#!/usr/bin/env python3
"""Select mitochondrial-read candidates from read-to-mitogenome PAF alignments."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--paf", required=True, type=Path)
    parser.add_argument("--read-ids", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--minimum-block-length", type=int, default=300)
    parser.add_argument("--minimum-identity", type=float, required=True)
    parser.add_argument("--minimum-union-query-bases", type=int, default=2000)
    return parser.parse_args()


def union_length(intervals: list[tuple[int, int]]) -> int:
    total = 0
    end = -1
    for start, stop in sorted(intervals):
        if stop <= end:
            continue
        total += stop - max(start, end)
        end = stop
    return total


def main() -> None:
    args = parse_args()
    intervals = defaultdict(list)
    read_lengths = {}
    with args.paf.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"Malformed PAF line {line_number} in {args.paf}")
            query_length = int(fields[1])
            query_start, query_end = int(fields[2]), int(fields[3])
            matches, block_length = int(fields[9]), int(fields[10])
            identity = matches / block_length if block_length else 0.0
            read_lengths[fields[0]] = query_length
            query_span = query_end - query_start
            if query_span >= args.minimum_block_length and identity >= args.minimum_identity:
                intervals[fields[0]].append((query_start, query_end))

    retained = []
    for read_id, read_intervals in intervals.items():
        covered = union_length(read_intervals)
        if covered >= args.minimum_union_query_bases:
            retained.append((read_id, read_lengths[read_id], covered, covered / read_lengths[read_id]))
    retained.sort()

    args.read_ids.parent.mkdir(parents=True, exist_ok=True)
    with args.read_ids.open("w", encoding="utf-8") as handle:
        for read_id, *_ in retained:
            handle.write(read_id + "\n")
    with args.summary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["read_id", "read_length", "mitochondrial_query_covered_bp", "covered_fraction"])
        for read_id, read_length, covered, fraction in retained:
            writer.writerow([read_id, read_length, covered, f"{fraction:.10f}"])
    print(f"Retained {len(retained)} mitochondrial-read candidates")


if __name__ == "__main__":
    main()
