import pytest
from src.mutation_caller import MutationCaller


def test_call_mutations_missense():
    wt_dna = "ATGCGACCCTCC"   # ATG (M), CGA (R), CCC (P), TCC (S)
    mut_dna = "ATGTGACCCTCC"  # ATG (M), TGA (*), CCC (P), TCC (S) -> Nonsense!

    mut_df = MutationCaller.call_mutations(wt_dna, mut_dna)

    assert len(mut_df) == 1
    assert mut_df.iloc[0]["position_nt"] == 4
    assert mut_df.iloc[0]["ref_nt"] == "C"
    assert mut_df.iloc[0]["mut_nt"] == "T"
    assert mut_df.iloc[0]["substitution_type"] == "Transition"
    assert "Nonsense" in mut_df.iloc[0]["functional_impact"]
