import pytest
from src.sequence_io import SequenceIOHandler


def test_dna_validation():
    assert SequenceIOHandler.validate_dna_sequence("ATGCGATCGATCG") is True
    assert SequenceIOHandler.validate_dna_sequence("ATGCX123") is False


def test_protein_validation():
    assert SequenceIOHandler.validate_protein_sequence("MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQ") is True
    assert SequenceIOHandler.validate_protein_sequence("MTEYK1234") is False
