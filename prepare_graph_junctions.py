#!/usr/bin/env python3
"""Create oriented junction and two-way branch definitions from a GFA graph."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gfa", required=True, type=Path)
    parser.add_argument("--junctions-output", required=True, type=Path)
    parser.add_argument("--branches-output", required=True, type=Path)
    parser.add_argument(
        "--id-map",
        type=Path,
        help="Optional TSV with original_analysis_id and new_junction_id columns",
    )
    return parser.parse_args()


def flip(state: str) -> str:
    return state[:-1] + ("-" if state[-1] == "+" else "+")


def read_id_map(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"original_analysis_id", "new_junction_id"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"{path} must contain: {', '.join(sorted(required))}")
    return {row["original_analysis_id"]: row["new_junction_id"] for row in rows}


def main() -> None:
    args = parse_args()
    id_map = read_id_map(args.id_map)
    links = []
    segment_names = set()
    with args.gfa.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if not fields or not fields[0]:
                continue
            if fields[0] == "S" and len(fields) >= 3:
                segment_names.add(fields[1])
            elif fields[0] == "L" and len(fields) >= 6:
                links.append(fields[1:6])
    if not links:
        raise SystemExit(f"No GFA L records found in {args.gfa}")

    records = []
    endpoint_junctions = defaultdict(set)
    for index, (from_name, from_orient, to_name, to_orient, overlap) in enumerate(links, 1):
        original_id = f"J{index:02d}"
        junction_id = id_map.get(original_id, original_id)
        from_state = from_name + from_orient
        to_state = to_name + to_orient
        records.append((junction_id, original_id, from_state, to_state, overlap))
        endpoint_junctions[from_state].add(junction_id)
        endpoint_junctions[flip(to_state)].add(junction_id)

    args.junctions_output.parent.mkdir(parents=True, exist_ok=True)
    with args.junctions_output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["junction_id", "original_analysis_id", "from_state", "to_state", "overlap"])
        writer.writerows(records)

    branch_sets = {}
    for endpoint, junction_ids in endpoint_junctions.items():
        if len(junction_ids) == 2:
            key = tuple(sorted(junction_ids))
            branch_sets.setdefault(key, endpoint)
    with args.branches_output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["branch_id", "oriented_endpoint", "junction_1", "junction_2"])
        for index, (junction_ids, endpoint) in enumerate(sorted(branch_sets.items()), 1):
            writer.writerow([f"B{index:02d}", endpoint, *junction_ids])

    print(f"Wrote {len(records)} junctions to {args.junctions_output}")
    print(f"Wrote {len(branch_sets)} two-way branches to {args.branches_output}")


if __name__ == "__main__":
    main()
