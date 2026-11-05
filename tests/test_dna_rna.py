import pytest
from src.dna_rna_analysis import calculate_gc_content, transcribe_and_translate, compute_sliding_window_gc


def test_calculate_gc_content():
    # Sequence with equal bases: 50%
    seq = "ATGC"
    assert calculate_gc_content(seq) == 50.0

    # All GC: 100%
    seq_gc = "GGCC"
    assert calculate_gc_content(seq_gc) == 100.0


def test_transcribe_and_translate():
    dna = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
    mrna, protein = transcribe_and_translate(dna)

    assert "U" in mrna
    assert "T" not in mrna
    assert protein.startswith("M")  # Methionine from ATG
