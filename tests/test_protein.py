import pytest
from src.protein_analysis import analyze_protein_properties


def test_analyze_protein_properties():
    # Poly-alanine peptide
    peptide = "AAAAAAAAAA"
    props = analyze_protein_properties(peptide)

    assert props["length_aa"] == 10
    assert props["molecular_weight_da"] > 700.0
    assert props["isoelectric_point_pi"] > 0.0
    assert props["gravy_hydropathy"] > 0.0  # Alanine is hydrophobic
