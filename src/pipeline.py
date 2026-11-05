import argparse
from pathlib import Path
from typing import Dict, Any, Tuple
import pandas as pd

from .sequence_io import SequenceIOHandler
from .dna_rna_analysis import calculate_gc_content, compute_sliding_window_gc, transcribe_and_translate
from .protein_analysis import analyze_protein_properties
from .alignment import SequenceAligner
from .mutation_caller import MutationCaller
from .phylogeny import PhylogeneticTreeBuilder
from .visualize import plot_gc_profile, plot_amino_acid_composition, plot_phylogenetic_tree


def run_bioinformatics_pipeline(
    variants_fasta: str | Path = "data/egfr_variants.fasta",
    homologs_fasta: str | Path = "data/tp53_homologs.fasta",
    output_dir: str | Path = "results",
) -> Dict[str, Any]:
    """Orchestrates comprehensive sequence analysis, variant identification, and phylogenetic modeling."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    # 1. Parse EGFR sequence variants
    egfr_records = SequenceIOHandler.read_fasta(variants_fasta)
    wt_record = egfr_records[0]
    t790m_record = egfr_records[1]
    l858r_record = egfr_records[2]

    # 2. DNA and Protein Characterization
    gc_wt = calculate_gc_content(wt_record.seq)
    mrna_wt, protein_wt = transcribe_and_translate(wt_record.seq)
    protein_props = analyze_protein_properties(protein_wt)

    # 3. Sliding Window GC & Skew Plot
    pos, gc_vals, skew_vals = compute_sliding_window_gc(wt_record.seq, window_size=40, step=5)
    plot_gc_profile(pos, gc_vals, skew_vals, fig_dir / "gc_skew_profile.png")

    # 4. Amino Acid Frequency Plot
    plot_amino_acid_composition(protein_props["aa_composition"], fig_dir / "amino_acid_profile.png")

    # 5. Pairwise Alignments
    aligner_nt = SequenceAligner(mode="global", is_protein=False)
    align_t790m = aligner_nt.align(str(wt_record.seq), str(t790m_record.seq))
    align_l858r = aligner_nt.align(str(wt_record.seq), str(l858r_record.seq))

    # 6. Mutation Calling (Clinically actionable oncogenic variants)
    mutations_t790m = MutationCaller.call_mutations(str(wt_record.seq), str(t790m_record.seq))
    mutations_l858r = MutationCaller.call_mutations(str(wt_record.seq), str(l858r_record.seq))
    all_mutations = pd.concat([
        mutations_t790m.assign(variant="EGFR_T790M"),
        mutations_l858r.assign(variant="EGFR_L858R"),
    ], ignore_index=True)
    all_mutations.to_csv(out_dir / "identified_mutations.csv", index=False)

    # 7. Homolog Phylogeny Construction
    tp53_records = SequenceIOHandler.read_fasta(homologs_fasta)
    tree_builder = PhylogeneticTreeBuilder(method="nj")
    tree, newick_str = tree_builder.build_tree_from_records(tp53_records)

    with open(out_dir / "tp53_phylogenetic_tree.nwk", "w", encoding="utf-8") as f:
        f.write(newick_str + "\n")

    plot_phylogenetic_tree(tree, fig_dir / "phylogenetic_tree.png")

    summary_report = {
        "wt_length_bp": len(wt_record.seq),
        "wt_gc_percent": gc_wt,
        "wt_protein_mw_da": protein_props["molecular_weight_da"],
        "wt_protein_pi": protein_props["isoelectric_point_pi"],
        "wt_gravy": protein_props["gravy_hydropathy"],
        "t790m_alignment_score": align_t790m["score"],
        "t790m_identity_pct": align_t790m["identity_percent"],
        "l858r_identity_pct": align_l858r["identity_percent"],
        "detected_mutations_count": len(all_mutations),
        "newick_tree": newick_str,
    }
    pd.DataFrame([summary_report]).to_csv(out_dir / "bioinformatics_summary.csv", index=False)
    return summary_report


def main():
    parser = argparse.ArgumentParser(description="Computational Bioinformatics Pipeline")
    parser.add_argument("--variants", "-v", type=str, default="data/egfr_variants.fasta", help="Variants FASTA")
    parser.add_argument("--homologs", "-m", type=str, default="data/tp53_homologs.fasta", help="Homologs FASTA")
    parser.add_argument("--output", "-o", type=str, default="results", help="Output directory")
    args = parser.parse_args()

    results = run_bioinformatics_pipeline(args.variants, args.homologs, output_dir=args.output)
    print("Bioinformatics Pipeline Execution Complete! Summary:")
    for k, v in results.items():
        if k != "newick_tree":
            print(f"  {k}: {v}")
    print(f"\nPhylogenetic Tree (Newick Format):\n  {results['newick_tree']}")


if __name__ == "__main__":
    main()
