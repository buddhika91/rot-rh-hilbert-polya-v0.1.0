from pathlib import Path
from rot_rh_hp.provenance import verify_sources

def test_recovered_research_sources_are_byte_integral():
    root = Path(__file__).resolve().parents[1]
    results = verify_sources(root)
    assert results
    assert all(results.values()), results
