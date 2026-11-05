import pytest
from src.alignment import SequenceAligner


def test_identical_sequence_alignment():
    aligner = SequenceAligner(mode="global", is_protein=False)
    seq = "ATGCGATCGATCG"
    res = aligner.align(seq, seq)

    assert res["identity_percent"] == 100.0
    assert res["matches"] == len(seq)


def test_protein_alignment_blosum62():
    aligner = SequenceAligner(mode="global", is_protein=True)
    seq1 = "HEAGAWGHEE"
    seq2 = "PAWHEAE"
    res = aligner.align(seq1, seq2)

    assert res["score"] > 0
    assert 0.0 < res["identity_percent"] <= 100.0
