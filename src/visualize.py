from pathlib import Path
from typing import List, Dict
import matplotlib.pyplot as plt
import numpy as np
from Bio import Phylo


def plot_gc_profile(
    positions: List[int],
    gc_values: List[float],
    skew_values: List[float],
    output_path: str | Path,
    title: str = "EGFR Coding Sequence: GC Content and Skew Analysis",
) -> Path:
    """Renders dual-axis plot of sliding-window GC percentage and GC skew."""
    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax1 = plt.subplots(figsize=(9, 5), dpi=300)

    # Primary axis: GC content
    color_gc = "#1565C0"
    ax1.set_xlabel("Nucleotide Sequence Position (bp)", fontsize=11, fontweight="bold")
    ax1.set_ylabel("GC Content (%)", color=color_gc, fontsize=11, fontweight="bold")
    line1 = ax1.plot(positions, gc_values, color=color_gc, linewidth=2.0, label="GC Content (%)")
    ax1.tick_params(axis="y", labelcolor=color_gc)
    ax1.axhline(50.0, color="gray", linestyle="--", linewidth=0.8, alpha=0.7)

    # Secondary axis: GC skew
    ax2 = ax1.twinx()
    color_skew = "#C2185B"
    ax2.set_ylabel("GC Skew (G-C)/(G+C)", color=color_skew, fontsize=11, fontweight="bold")
    line2 = ax2.plot(positions, skew_values, color=color_skew, linewidth=1.8, linestyle=":", label="GC Skew")
    ax2.tick_params(axis="y", labelcolor=color_skew)
    ax2.axhline(0.0, color="black", linestyle="-", linewidth=0.6, alpha=0.5)

    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc="upper right", fontsize=9)

    plt.title(title, fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()
    return out_path


def plot_amino_acid_composition(
    aa_comp: Dict[str, float],
    output_path: str | Path,
    title: str = "Protein Residue Abundance Spectrum",
) -> Path:
    """Plots amino acid composition percentages with chemical class color grouping."""
    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    labels = list(aa_comp.keys())
    values = list(aa_comp.values())

    # Residue classification
    pos_charged = {"K", "R", "H"}
    neg_charged = {"D", "E"}
    aromatic = {"F", "Y", "W"}

    colors = []
    for aa in labels:
        if aa in pos_charged:
            colors.append("#1E88E5")  # Blue for basic
        elif aa in neg_charged:
            colors.append("#E53935")  # Red for acidic
        elif aa in aromatic:
            colors.append("#8E24AA")  # Purple for aromatic
        else:
            colors.append("#43A047")  # Green for neutral/aliphatic

    fig, ax = plt.subplots(figsize=(10, 5), dpi=300)
    bars = ax.bar(labels, values, color=colors, edgecolor="black", alpha=0.85)

    ax.set_xlabel("Amino Acid Residue", fontsize=11, fontweight="bold")
    ax.set_ylabel("Frequency (%)", fontsize=11, fontweight="bold")
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.grid(axis="y", linestyle=":", alpha=0.6)

    for bar in bars:
        h = bar.get_height()
        if h > 0:
            ax.text(bar.get_x() + bar.get_width()/2, h + 0.15, f"{h:.1f}", ha="center", va="bottom", fontsize=7.5)

    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()
    return out_path


def plot_phylogenetic_tree(
    tree: Phylo.BaseTree.Tree,
    output_path: str | Path,
    title: str = "Phylogenetic Evolutionary Tree (Neighbor-Joining)",
) -> Path:
    """Renders phylogenetic tree dendrogram to image file."""
    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(9, 6), dpi=300)
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, do_show=False, axes=ax)
    ax.set_title(title, fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()
    return out_path
