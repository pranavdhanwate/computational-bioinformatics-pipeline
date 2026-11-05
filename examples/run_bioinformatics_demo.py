import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

from src.pipeline import run_bioinformatics_pipeline

if __name__ == "__main__":
    variants_fa = project_root / "data" / "egfr_variants.fasta"
    homologs_fa = project_root / "data" / "tp53_homologs.fasta"
    out_dir = project_root / "results"

    print("Running Computational Bioinformatics Sequence Analysis Demo...")
    summary = run_bioinformatics_pipeline(variants_fa, homologs_fa, output_dir=out_dir)

    print("\nPipeline Summary Results:")
    for k, v in summary.items():
        if k != "newick_tree":
            print(f"  {k}: {v}")
    print(f"\nPhylogenetic Tree (Newick): {summary['newick_tree']}")
