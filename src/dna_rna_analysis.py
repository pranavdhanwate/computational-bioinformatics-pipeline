from typing import List, Dict, Tuple, Any
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction


def calculate_gc_content(sequence: str | Seq) -> float:
    """Calculates GC percentage of a DNA/RNA nucleotide sequence."""
    seq_str = str(sequence).upper()
    total_len = len(seq_str)
    if total_len == 0:
        return 0.0
    gc_count = seq_str.count("G") + seq_str.count("C")
    return round((gc_count / total_len) * 100.0, 2)


def compute_sliding_window_gc(sequence: str | Seq, window_size: int = 50, step: int = 10) -> Tuple[List[int], List[float], List[float]]:
    """Calculates sliding-window GC content (%) and GC skew: (G - C) / (G + C)."""
    seq_str = str(sequence).upper()
    positions = []
    gc_values = []
    skew_values = []

    for i in range(0, len(seq_str) - window_size + 1, step):
        window = seq_str[i : i + window_size]
        g = window.count("G")
        c = window.count("C")
        total_gc = g + c

        gc_pct = (total_gc / window_size) * 100.0
        skew = (g - c) / total_gc if total_gc > 0 else 0.0

        positions.append(i + window_size // 2)
        gc_values.append(round(gc_pct, 2))
        skew_values.append(round(skew, 4))

    return positions, gc_values, skew_values


def transcribe_and_translate(dna_sequence: str | Seq, table: int = 1) -> Tuple[str, str]:
    """Transcribes coding DNA to mRNA and translates mRNA to peptide using standard genetic code."""
    seq_obj = Seq(str(dna_sequence).upper())
    mrna = seq_obj.transcribe()
    # Trim to multiple of 3 for translation
    remainder = len(mrna) % 3
    translatable_mrna = mrna[:-remainder] if remainder > 0 else mrna
    protein = translatable_mrna.translate(table=table, to_stop=False)
    return str(mrna), str(protein)


def find_open_reading_frames(dna_sequence: str | Seq, min_aa_length: int = 30) -> List[Dict[str, Any]]:
    """Identifies Open Reading Frames (ORFs) from start codon (ATG) to stop (TAA, TAG, TGA) across all 6 frames."""
    from typing import Any
    seq_obj = Seq(str(dna_sequence).upper())
    orfs = []

    for strand, nuc_seq in [(+1, seq_obj), (-1, seq_obj.reverse_complement())]:
        for frame in range(3):
            sub_seq = nuc_seq[frame:]
            trans = sub_seq.translate(to_stop=False)
            aa_str = str(trans)
            # Find ORFs starting with 'M' and ending with '*'
            start_idx = 0
            while True:
                m_pos = aa_str.find("M", start_idx)
                if m_pos == -1:
                    break
                stop_pos = aa_str.find("*", m_pos)
                if stop_pos == -1:
                    break
                protein_len = stop_pos - m_pos
                if protein_len >= min_aa_length:
                    peptide = aa_str[m_pos:stop_pos]
                    start_nt = frame + (m_pos * 3) if strand == +1 else len(seq_obj) - (frame + stop_pos * 3)
                    end_nt = frame + ((stop_pos + 1) * 3) if strand == +1 else len(seq_obj) - (frame + m_pos * 3)
                    orfs.append({
                        "strand": "+" if strand == +1 else "-",
                        "frame": frame + 1,
                        "start_nt": start_nt,
                        "end_nt": end_nt,
                        "length_aa": protein_len,
                        "protein_sequence": peptide,
                    })
                start_idx = m_pos + 1

    return sorted(orfs, key=lambda x: x["length_aa"], reverse=True)
