import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "call_mitogenome_variants.py"
SPEC = importlib.util.spec_from_file_location("variant_caller", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_snv_and_indels_are_coordinate_aware():
    calls = MODULE.variants_from_alignment(
        query_aligned="ATGCA-T",
        reference_aligned="AC-CAGT",
        query_start=101,
        query_end=106,
        reference_start=201,
        strand="+",
        block_id="block_1",
    )
    assert [call["variant_type"] for call in calls] == ["SNV", "insertion", "deletion"]
    assert calls[0]["reference_start"] == 202
    assert calls[0]["query_start"] == 102
    assert calls[1]["alt"] == "G"
    assert calls[2]["ref"] == "G"


def test_reverse_query_coordinates_decrease():
    calls = MODULE.variants_from_alignment(
        query_aligned="AT",
        reference_aligned="AG",
        query_start=10,
        query_end=11,
        reference_start=20,
        strand="-",
        block_id="block_1",
    )
    assert calls[0]["query_start"] == 10
    assert calls[0]["reference_start"] == 21
