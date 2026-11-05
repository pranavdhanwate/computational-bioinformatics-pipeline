from typing import List, Dict, Any
import pandas as pd
from Bio.Seq import Seq


class MutationCaller:
    """Identifies nucleotide polymorphisms and translates functional impact on protein sequence."""

    PURINES = {"A", "G"}
    PYRIMIDINES = {"C", "T"}

    @classmethod
    def classify_nucleotide_substitution(cls, ref_nt: str, mut_nt: str) -> str:
        ref = ref_nt.upper()
        mut = mut_nt.upper()
        if (ref in cls.PURINES and mut in cls.PURINES) or (ref in cls.PYRIMIDINES and mut in cls.PYRIMIDINES):
            return "Transition"
        return "Transversion"

    @classmethod
    def call_mutations(cls, wild_type_dna: str, mutant_dna: str) -> pd.DataFrame:
        wt = str(wild_type_dna).upper().strip()
        mut = str(mutant_dna).upper().strip()
        min_len = min(len(wt), len(mut))

        mutations = []
        for i in range(min_len):
            if wt[i] != mut[i]:
                pos_1based = i + 1
                ref_nt = wt[i]
                mut_nt = mut[i]
                sub_type = cls.classify_nucleotide_substitution(ref_nt, mut_nt)

                # Identify codon position
                codon_idx = i // 3
                codon_start = codon_idx * 3
                codon_end = codon_start + 3

                wt_codon = wt[codon_start:codon_end] if codon_end <= len(wt) else ""
                mut_codon = mut[codon_start:codon_end] if codon_end <= len(mut) else ""

                wt_aa = str(Seq(wt_codon).translate()) if len(wt_codon) == 3 else "?"
                mut_aa = str(Seq(mut_codon).translate()) if len(mut_codon) == 3 else "?"

                if wt_aa == mut_aa:
                    impact = "Synonymous (Silent)"
                elif mut_aa == "*":
                    impact = "Nonsense (Premature Stop)"
                else:
                    impact = f"Missense ({wt_aa}{codon_idx + 1}{mut_aa})"

                mutations.append({
                    "position_nt": pos_1based,
                    "codon_number": codon_idx + 1,
                    "ref_nt": ref_nt,
                    "mut_nt": mut_nt,
                    "substitution_type": sub_type,
                    "wt_codon": wt_codon,
                    "mut_codon": mut_codon,
                    "wt_amino_acid": wt_aa,
                    "mut_amino_acid": mut_aa,
                    "functional_impact": impact,
                })

        return pd.DataFrame(mutations)
