#!/usr/bin/env python3
"""Summarize assembly-graph junction support from read-level traversal events.

The input contains one row per traversal event. A read can therefore occur more
than once at the same junction. The script reports both unique supporting reads
and total traversal events, and computes two-way branch composition with
whole-read cluster bootstrap confidence intervals.
"""

from __future__ import annotations

import argparse
import csv
import random
from collections import Counter, defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--events", required=True, type=Path, help="Headered event-level TSV")
    parser.add_argument("--branches", required=True, type=Path, help="Headered branch-definition TSV")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--minimum-graph-query-coverage", type=float, default=0.95)
    parser.add_argument(
        "--exclude-provenance", default="nuclear_boundary_NUMT",
        help="Comma-separated provenance classes to exclude",
    )
    parser.add_argument("--bootstrap-replicates", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260812)
    return parser.parse_args()


def read_tsv(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def require_columns(path: Path, rows, required):
    if not rows:
        raise ValueError(f"No data rows found in {path}")
    missing = set(required) - set(rows[0])
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")


def percentile(values, probability):
    values = sorted(values)
    if not values:
        return float("nan")
    index = (len(values) - 1) * probability
    lower = int(index)
    upper = min(lower + 1, len(values) - 1)
    fraction = index - lower
    return values[lower] * (1 - fraction) + values[upper] * fraction


def bootstrap_fraction(read_events, junction_1, junction_2, replicates, rng):
    read_ids = list(read_events)
    estimates = []
    for _ in range(replicates):
        counts = Counter()
        for _ in read_ids:
            sampled_read = rng.choice(read_ids)
            counts.update(read_events[sampled_read])
        denominator = counts[junction_1] + counts[junction_2]
        if denominator:
            estimates.append(counts[junction_1] / denominator)
    return percentile(estimates, 0.025), percentile(estimates, 0.975)


def main() -> None:
    args = parse_args()
    events = read_tsv(args.events)
    branches = read_tsv(args.branches)
    require_columns(
        args.events, events,
        {"sample", "read_id", "junction_id", "graph_query_coverage", "provenance_class"},
    )
    require_columns(
        args.branches, branches,
        {"branch_id", "oriented_endpoint", "junction_1", "junction_2"},
    )
    excluded = {value for value in args.exclude_provenance.split(",") if value}
    retained = [
        row for row in events
        if float(row["graph_query_coverage"]) >= args.minimum_graph_query_coverage
        and row["provenance_class"] not in excluded
    ]
    if not retained:
        raise SystemExit("No traversal events remain after filtering")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    by_sample_read = defaultdict(lambda: defaultdict(list))
    for row in retained:
        by_sample_read[row["sample"]][row["read_id"]].append(row["junction_id"])

    junctions = sorted({row["junction_id"] for row in retained})
    junction_path = args.output_dir / "junction_support.tsv"
    with junction_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "junction_id", "distinct_reads", "traversal_events"])
        for sample in sorted(by_sample_read):
            for junction in junctions:
                supporting_reads = sum(
                    junction in read_events for read_events in by_sample_read[sample].values()
                )
                traversal_events = sum(
                    read_events.count(junction) for read_events in by_sample_read[sample].values()
                )
                writer.writerow([sample, junction, supporting_reads, traversal_events])

    rng = random.Random(args.seed)
    branch_path = args.output_dir / "branch_frequencies.tsv"
    with branch_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "sample", "branch_id", "oriented_endpoint", "junction_1", "junction_2",
            "junction_1_distinct_reads", "junction_2_distinct_reads", "reads_supporting_both",
            "junction_1_distinct_read_fraction", "junction_2_distinct_read_fraction",
            "junction_1_events", "junction_2_events", "total_events",
            "junction_1_event_fraction", "junction_2_event_fraction",
            "junction_1_event_fraction_read_cluster_bootstrap95_low",
            "junction_1_event_fraction_read_cluster_bootstrap95_high",
        ])
        for sample in sorted(by_sample_read):
            read_events = by_sample_read[sample]
            for branch in branches:
                j1, j2 = branch["junction_1"], branch["junction_2"]
                j1_reads = sum(j1 in values for values in read_events.values())
                j2_reads = sum(j2 in values for values in read_events.values())
                both = sum(j1 in values and j2 in values for values in read_events.values())
                distinct_total = j1_reads + j2_reads
                distinct_fraction_1 = j1_reads / distinct_total if distinct_total else float("nan")
                distinct_fraction_2 = j2_reads / distinct_total if distinct_total else float("nan")
                j1_events = sum(values.count(j1) for values in read_events.values())
                j2_events = sum(values.count(j2) for values in read_events.values())
                total = j1_events + j2_events
                fraction_1 = j1_events / total if total else float("nan")
                fraction_2 = j2_events / total if total else float("nan")
                low, high = bootstrap_fraction(
                    read_events, j1, j2, args.bootstrap_replicates, rng
                ) if total else (float("nan"), float("nan"))
                writer.writerow([
                    sample, branch["branch_id"], branch["oriented_endpoint"], j1, j2,
                    j1_reads, j2_reads, both,
                    f"{distinct_fraction_1:.10f}", f"{distinct_fraction_2:.10f}",
                    j1_events, j2_events, total,
                    f"{fraction_1:.10f}", f"{fraction_2:.10f}",
                    f"{low:.10f}", f"{high:.10f}",
                ])

    print(f"Retained {len(retained)} traversal events")
    print(f"Wrote {junction_path}")
    print(f"Wrote {branch_path}")


if __name__ == "__main__":
    main()
