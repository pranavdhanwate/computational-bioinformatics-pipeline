import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from src.phylogeny import PhylogeneticTreeBuilder


def test_build_phylogenetic_tree():
    rec1 = SeqRecord(Seq("ATGCGATCGATCGATC"), id="SpeciesA")
    rec2 = SeqRecord(Seq("ATGCGATCGATCGATT"), id="SpeciesB")
    rec3 = SeqRecord(Seq("TTGCGATCGATCGATT"), id="SpeciesC")

    builder = PhylogeneticTreeBuilder(method="nj")
    tree, newick = builder.build_tree_from_records([rec1, rec2, rec3])

    assert tree is not None
    assert len(newick) > 0
    assert "SpeciesA" in newick
    assert "SpeciesB" in newick
