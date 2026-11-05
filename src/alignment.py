from typing import Dict, Any, Tuple
from Bio.Align import PairwiseAligner, substitution_matrices


class SequenceAligner:
    """Performs global (Needleman-Wunsch) and local (Smith-Waterman) sequence alignments."""

    def __init__(self, mode: str = "global", is_protein: bool = False):
        self.aligner = PairwiseAligner()
        self.aligner.mode = mode
        self.is_protein = is_protein

        if is_protein:
            self.aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
            self.aligner.open_gap_score = -10.0
            self.aligner.extend_gap_score = -0.5
        else:
            self.aligner.match_score = 2.0
            self.aligner.mismatch_score = -1.0
            self.aligner.open_gap_score = -5.0
            self.aligner.extend_gap_score = -1.0

    def align(self, seq_a: str, seq_b: str) -> Dict[str, Any]:
        """Aligns two biological sequences and reports score and identity percentage."""
        str_a = str(seq_a).strip().upper()
        str_b = str(seq_b).strip().upper()

        alignments = list(self.aligner.align(str_a, str_b))
        if not alignments:
            return {"score": 0.0, "identity_percent": 0.0, "alignment_str": "No alignment"}

        best_aln = alignments[0]
        score = float(best_aln.score)

        # Compute exact identity count
        aligned_pairs = list(zip(best_aln[0], best_aln[1]))
        matches = sum(1 for a, b in aligned_pairs if a == b and a != "-")
        total_len = max(len(best_aln[0]), 1)
        identity_pct = (matches / total_len) * 100.0

        return {
            "mode": self.aligner.mode,
            "score": round(score, 2),
            "identity_percent": round(identity_pct, 2),
            "matches": matches,
            "alignment_length": total_len,
            "aligned_seq_a": str(best_aln[0]),
            "aligned_seq_b": str(best_aln[1]),
            "alignment_representation": str(best_aln),
        }
