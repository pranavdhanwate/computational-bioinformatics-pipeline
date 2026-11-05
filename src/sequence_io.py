from pathlib import Path
from typing import List, Dict, Any, Generator
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class SequenceIOHandler:
    """Handles parsing, validation, and serialization of biological FASTA files."""

    @staticmethod
    def read_fasta(filepath: str | Path) -> List[SeqRecord]:
        path = Path(filepath)
        if not path.exists():
            raise FileNotFoundError(f"FASTA file not found: {path}")

        records = list(SeqIO.parse(str(path), "fasta"))
        if not records:
            raise ValueError(f"No valid FASTA records found in {path}")
        return records

    @staticmethod
    def validate_dna_sequence(sequence: str) -> bool:
        valid_bases = set("ATCGNatcgn")
        return set(str(sequence)).issubset(valid_bases)

    @staticmethod
    def validate_protein_sequence(sequence: str) -> bool:
        valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY*acdefghiklmnpqrstvwy*")
        return set(str(sequence)).issubset(valid_amino_acids)

    @staticmethod
    def write_fasta(records: List[SeqRecord], output_path: str | Path) -> Path:
        out_path = Path(output_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        SeqIO.write(records, str(out_path), "fasta")
        return out_path
