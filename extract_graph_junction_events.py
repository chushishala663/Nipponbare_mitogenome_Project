#!/usr/bin/env python3
"""Convert whole-read GAF paths and nuclear PAF alignments to junction events."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Alignment:
    query_start: int
    query_end: int
    target: str
    target_start: int
    target_end: int
    identity: float
    path: tuple[str, ...] = ()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--gaf", required=True, type=Path)
    parser.add_argument("--nuclear-paf", required=True, type=Path)
    parser.add_argument("--numt-bed", required=True, type=Path)
    parser.add_argument("--junctions", required=True, type=Path)
    parser.add_argument("--events-output", required=True, type=Path)
    parser.add_argument("--provenance-output", required=True, type=Path)
    parser.add_argument("--minimum-graph-identity", type=float, required=True)
    parser.add_argument("--minimum-nuclear-identity", type=float, default=0.70)
    parser.add_argument("--maximum-alignment-overlap", type=int, default=500)
    parser.add_argument("--maximum-split-query-gap", type=int, default=2000)
    parser.add_argument("--maximum-split-target-gap", type=int, default=10000)
    parser.add_argument("--nuclear-flank", type=int, default=5000)
    parser.add_argument("--internal-numt-query-fraction", type=float, default=0.80)
    return parser.parse_args()


def flip(state: str) -> str:
    return state[:-1] + ("-" if state[-1] == "+" else "+")


def reverse_path(path: tuple[str, ...]) -> tuple[str, ...]:
    return tuple(flip(state) for state in reversed(path))


def union_length(intervals: list[tuple[int, int]]) -> int:
    total = 0
    end = -1
    for start, stop in sorted(intervals):
        if stop <= end:
            continue
        total += stop - max(start, end)
        end = stop
    return total


def tag_value(fields: list[str], prefix: str) -> str:
    for field in fields[12:]:
        if field.startswith(prefix):
            return field[len(prefix):]
    return ""


def parse_gaf(path: Path, minimum_identity: float) -> tuple[dict[str, int], dict[str, list[Alignment]]]:
    lengths = {}
    by_read = defaultdict(list)
    path_pattern = re.compile(r"([><])([^><]+)")
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"Malformed GAF line {line_number} in {path}")
            read_id = fields[0]
            lengths[read_id] = int(fields[1])
            matches, block_length = int(fields[9]), int(fields[10])
            identity_tag = tag_value(fields, "id:f:")
            identity = float(identity_tag) if identity_tag else (matches / block_length if block_length else 0.0)
            if identity < minimum_identity:
                continue
            states = tuple(name + ("+" if marker == ">" else "-") for marker, name in path_pattern.findall(fields[5]))
            if not states:
                continue
            if fields[4] == "-":
                states = reverse_path(states)
            by_read[read_id].append(
                Alignment(int(fields[2]), int(fields[3]), fields[5], int(fields[7]), int(fields[8]), identity, states)
            )
    return lengths, by_read


def select_query_chain(alignments: list[Alignment], maximum_overlap: int) -> list[Alignment]:
    selected = []
    for alignment in sorted(
        alignments,
        key=lambda item: ((item.query_end - item.query_start) * item.identity, item.identity),
        reverse=True,
    ):
        overlap = max(
            (min(alignment.query_end, kept.query_end) - max(alignment.query_start, kept.query_start) for kept in selected),
            default=0,
        )
        if overlap <= maximum_overlap:
            selected.append(alignment)
    return sorted(selected, key=lambda item: (item.query_start, item.query_end))


def read_junctions(path: Path) -> dict[tuple[str, str], str]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"junction_id", "from_state", "to_state"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"{path} must contain: {', '.join(sorted(required))}")
    arcs = {}
    for row in rows:
        source, target = row["from_state"], row["to_state"]
        arcs[(source, target)] = row["junction_id"]
        arcs[(flip(target), flip(source))] = row["junction_id"]
    return arcs


def parse_nuclear_paf(path: Path, minimum_identity: float) -> dict[str, list[Alignment]]:
    by_read = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"Malformed PAF line {line_number} in {path}")
            block_length = int(fields[10])
            identity = int(fields[9]) / block_length if block_length else 0.0
            if identity >= minimum_identity:
                by_read[fields[0]].append(
                    Alignment(int(fields[2]), int(fields[3]), fields[5], int(fields[7]), int(fields[8]), identity)
                )
    return by_read


def read_numt_bed(path: Path) -> list[tuple[str, int, int, str]]:
    loci = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"Malformed BED line {line_number} in {path}")
            loci.append((fields[0], int(fields[1]), int(fields[2]), fields[3] if len(fields) > 3 else f"NUMT{line_number:03d}"))
    return loci


def intervals_linked(a: Alignment, b: Alignment, max_query_gap: int, max_target_gap: int) -> bool:
    query_gap = max(a.query_start, b.query_start) - min(a.query_end, b.query_end)
    target_gap = max(a.target_start, b.target_start) - min(a.target_end, b.target_end)
    return a.target == b.target and query_gap <= max_query_gap and target_gap <= max_target_gap


def classify_provenance(
    alignments: list[Alignment],
    read_length: int,
    loci: list[tuple[str, int, int, str]],
    nuclear_flank: int,
    internal_fraction: float,
    max_query_gap: int,
    max_target_gap: int,
) -> tuple[str, str, float, int]:
    best_internal = (0.0, "")
    partial_locus = ""
    for chrom, start, end, locus_id in loci:
        relevant = [a for a in alignments if a.target == chrom and a.target_end > start and a.target_start < end]
        if not relevant:
            continue
        partial_locus = partial_locus or locus_id
        linked = list(relevant)
        for overlap_alignment in relevant:
            linked.extend(
                a for a in alignments
                if intervals_linked(overlap_alignment, a, max_query_gap, max_target_gap)
            )
        linked = list(set(linked))
        left_extension = max((start - a.target_start for a in linked if a.target_start < start), default=0)
        right_extension = max((a.target_end - end for a in linked if a.target_end > end), default=0)
        if max(left_extension, right_extension) >= nuclear_flank:
            covered = union_length([(a.query_start, a.query_end) for a in linked])
            return "nuclear_boundary_NUMT", locus_id, covered / read_length, max(left_extension, right_extension)
        covered = union_length([(a.query_start, a.query_end) for a in relevant])
        fraction = covered / read_length
        if fraction > best_internal[0]:
            best_internal = (fraction, locus_id)
    if best_internal[0] >= internal_fraction:
        return "internal_NUMT_ambiguous", best_internal[1], best_internal[0], 0
    if partial_locus:
        return "partial_NUMT_similarity", partial_locus, best_internal[0], 0
    return "other_retained", "", 0.0, 0


def main() -> None:
    args = parse_args()
    junctions = read_junctions(args.junctions)
    read_lengths, graph_alignments = parse_gaf(args.gaf, args.minimum_graph_identity)
    nuclear_alignments = parse_nuclear_paf(args.nuclear_paf, args.minimum_nuclear_identity)
    loci = read_numt_bed(args.numt_bed)

    args.events_output.parent.mkdir(parents=True, exist_ok=True)
    provenance_rows = []
    event_rows = []
    for read_id in sorted(graph_alignments):
        chain = select_query_chain(graph_alignments[read_id], args.maximum_alignment_overlap)
        read_length = read_lengths[read_id]
        graph_covered = union_length([(a.query_start, a.query_end) for a in chain])
        graph_coverage = graph_covered / read_length
        provenance, locus_id, nuclear_fraction, flank_extension = classify_provenance(
            nuclear_alignments.get(read_id, []), read_length, loci,
            args.nuclear_flank, args.internal_numt_query_fraction,
            args.maximum_split_query_gap, args.maximum_split_target_gap,
        )
        provenance_rows.append([
            args.sample, read_id, read_length, graph_covered, f"{graph_coverage:.10f}",
            provenance, locus_id, f"{nuclear_fraction:.10f}", flank_extension,
        ])
        event_index = 0
        for alignment_index, alignment in enumerate(chain):
            for source, target in zip(alignment.path, alignment.path[1:]):
                junction_id = junctions.get((source, target))
                if junction_id:
                    event_index += 1
                    event_rows.append([
                        args.sample, read_id, event_index, junction_id, source, target,
                        "within_gaf_path", f"{graph_coverage:.10f}", provenance,
                    ])
            if alignment_index + 1 >= len(chain):
                continue
            following = chain[alignment_index + 1]
            query_gap = following.query_start - alignment.query_end
            if query_gap <= args.maximum_split_query_gap and query_gap >= -args.maximum_alignment_overlap:
                source, target = alignment.path[-1], following.path[0]
                junction_id = junctions.get((source, target))
                if junction_id:
                    event_index += 1
                    event_rows.append([
                        args.sample, read_id, event_index, junction_id, source, target,
                        "linked_split_gaf", f"{graph_coverage:.10f}", provenance,
                    ])

    with args.provenance_output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "sample", "read_id", "read_length", "graph_query_covered_bp",
            "graph_query_coverage", "provenance_class", "numt_locus",
            "nuclear_numt_query_fraction", "nuclear_flank_extension_bp",
        ])
        writer.writerows(provenance_rows)
    with args.events_output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "sample", "read_id", "event_index", "junction_id", "from_state", "to_state",
            "event_source", "graph_query_coverage", "provenance_class",
        ])
        writer.writerows(event_rows)
    print(f"Wrote {len(provenance_rows)} read classifications to {args.provenance_output}")
    print(f"Wrote {len(event_rows)} traversal events to {args.events_output}")


if __name__ == "__main__":
    main()
