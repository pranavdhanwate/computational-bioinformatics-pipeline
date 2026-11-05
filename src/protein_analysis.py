from typing import Dict, Any
from Bio.SeqUtils.ProtParam import ProteinAnalysis


def analyze_protein_properties(protein_sequence: str) -> Dict[str, Any]:
    """Computes physicochemical and structural characteristics of a peptide sequence."""
    clean_seq = str(protein_sequence).replace("*", "").replace("X", "").strip()
    if not clean_seq:
        return {}

    analyst = ProteinAnalysis(clean_seq)
    mw = analyst.molecular_weight()
    pi = analyst.isoelectric_point()
    gravy = analyst.gravy()
    aromaticity = analyst.aromaticity()
    sec_struc = analyst.secondary_structure_fraction()  # (Helix, Turn, Sheet)
    aa_percent = analyst.amino_acids_percent

    # Basic vs Acidic residues
    basic_count = clean_seq.count("K") + clean_seq.count("R") + clean_seq.count("H")
    acidic_count = clean_seq.count("D") + clean_seq.count("E")

    return {
        "length_aa": len(clean_seq),
        "molecular_weight_da": round(float(mw), 2),
        "isoelectric_point_pi": round(float(pi), 2),
        "gravy_hydropathy": round(float(gravy), 3),
        "aromaticity": round(float(aromaticity), 3),
        "helix_fraction": round(float(sec_struc[0]), 3),
        "turn_fraction": round(float(sec_struc[1]), 3),
        "sheet_fraction": round(float(sec_struc[2]), 3),
        "basic_residues": int(basic_count),
        "acidic_residues": int(acidic_count),
        "aa_composition": {aa: round(pct * 100.0, 2) for aa, pct in aa_percent.items()},
    }
