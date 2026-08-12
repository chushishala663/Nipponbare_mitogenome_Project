import csv
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_distinct_reads_and_events_are_counted_separately():
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp)
        subprocess.run([
            sys.executable, str(ROOT / "scripts" / "summarize_graph_junctions.py"),
            "--events", str(ROOT / "examples" / "junction_events.example.tsv"),
            "--branches", str(ROOT / "examples" / "branches.example.tsv"),
            "--output-dir", str(out), "--bootstrap-replicates", "20",
        ], check=True)
        junctions = rows(out / "junction_support.tsv")
        hifi_j01 = next(row for row in junctions if row["sample"] == "HIFI" and row["junction_id"] == "J01")
        ont_j01 = next(row for row in junctions if row["sample"] == "ONT" and row["junction_id"] == "J01")
        assert hifi_j01["distinct_reads"] == "2"
        assert hifi_j01["traversal_events"] == "3"
        assert ont_j01["distinct_reads"] == "1"
        assert ont_j01["traversal_events"] == "2"
        branches = rows(out / "branch_frequencies.tsv")
        ont = next(row for row in branches if row["sample"] == "ONT")
        assert ont["junction_1_events"] == "2"
        assert ont["junction_2_events"] == "1"
        assert ont["total_events"] == "3"
