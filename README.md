# Computational Bioinformatics Pipeline

A biological sequence analysis, mutation identification, and molecular evolution framework built using Biopython for genomic and proteomic investigation of therapeutic targets and disease variants.

---

## 1. Project Title
**Computational Bioinformatics Pipeline: Sequence Analysis, Mutation Detection, and Phylogenetic Reconstruction**

## 2. Short Scientific Description
This project provides a biological sequence processing pipeline for DNA, RNA, and protein sequences. The system executes sequence validation, nucleotide composition and GC skew profiling, open reading frame (ORF) extraction, in-silico translation, ProtParam physicochemical profiling, pairwise global and local alignments using substitution matrices (BLOSUM62/PAM250), automated mutation detection with clinical nomenclature, and phylogenetic reconstruction via distance-based clustering algorithms (UPGMA and Neighbor-Joining).

## 3. Problem Statement
Translating raw next-generation sequencing (NGS) data into clinically actionable therapeutic hypotheses requires automated sequence analysis. Researchers frequently encounter non-standard fasta formatting, unverified IUPAC characters, fragmented alignment scripts, and manual variant curation bottlenecks that limit reproducible comparative genomic analysis across therapeutic targets.

## 4. Scientific Relevance
Accurate variant detection is essential for predicting acquired drug resistance in targeted oncology. This pipeline demonstrates variant discovery using the human **Epidermal Growth Factor Receptor (EGFR)** kinase sequence, detecting the classical gatekeeper mutations **T790M** (gefitinib/erlotinib resistance) and **L858R** (activating driver mutation). Evolutionary divergence is demonstrated through phylogenetic analysis of the **TP53** tumor suppressor gene across six mammalian and vertebrate homologues (human, chimpanzee, rhesus macaque, mouse, rat, and dog).

## 5. Key Capabilities
- **Biological Sequence Parsing:** Multi-FASTA file ingestion with strict IUPAC character verification and validation for DNA, RNA, and amino acid alphabets.
- **Nucleotide Architecture Profiling:** Base composition analysis, overall GC content calculation, and sliding-window GC skew ($[G-C]/[G+C]$) profiling to identify replication origins and transcription boundaries.
- **In-Silico Expression Pipeline:** DNA transcription into mRNA, reading frame detection, and codon-based translation into corresponding polypeptide sequences.
- **Protein Physicochemical Profiling:** Calculation of molecular weight, theoretical isoelectric point (pI), GRAVY hydrophobicity index, aromaticity index, instability index, and secondary structure fraction.
- **Pairwise Sequence Alignment:** Implementation of global Needleman-Wunsch and local Smith-Waterman algorithms supporting affine gap penalties and substitution matrices (BLOSUM62, PAM250).
- **Automated Mutation Calling:** Alignment-based variant identification reporting 1-based indexing, wild-type vs. mutant residues, and standard clinical nomenclature (e.g., p.Thr790Met).
- **Phylogenetic Reconstruction:** Calculation of identity distance matrices, UPGMA hierarchical clustering, and Saitou-Nei Neighbor-Joining (NJ) tree generation with Newick export.
- **Automated Visualisation:** Sliding-window GC skew graphs, mutation distribution plots, and rooted phylogenetic trees.

## 6. Technology Stack
- **Language:** Python 3.10+
- **Bioinformatics Engine:** Biopython (1.83+)
- **External Alignment Tools:** BLAST+ (optional external CLI), MAFFT (optional)
- **Data & Numerical Processing:** NumPy, Pandas
- **Visualization:** Matplotlib
- **Test Framework:** Pytest

## 7. Architecture
The architecture emphasizes modularity and biological type safety:
- src/sequence_io.py: Multi-FASTA ingestion and strict IUPAC biological character validation.
- src/dna_rna_analysis.py: Nucleotide composition, GC content, sliding-window GC skew, transcription, and translation.
- src/protein_analysis.py: ProtParam physicochemical properties (MW, pI, GRAVY, aromaticity, instability).
- src/alignment.py: Pairwise global/local alignment with BLOSUM62 scoring and gap penalties.
- src/mutation_caller.py: Sequence alignment variant detection with 1-based clinical nomenclature.
- src/phylogeny.py: Pairwise distance matrix computation, UPGMA and Neighbor-Joining tree construction.
- src/visualize.py: GC skew plots, mutation distribution diagrams, and phylogenetic tree rendering.
- src/pipeline.py: Comprehensive automated workflow execution.

## 8. Workflow Diagram
`	ext
data/egfr_variants.fasta ────────┐
data/tp53_homologs.fasta ───────┐│
                                ▼▼
                     [1. Sequence IO & Validation]
                     (Multi-FASTA parsing, IUPAC validation)
                                │
                                ▼
                     [2. Nucleotide & Protein Analysis]
                     (GC content, skew, translation, ProtParam)
                                │
                                ▼
                     [3. Pairwise / Multiple Alignment]
                     (Needleman-Wunsch, BLOSUM62, gap penalties)
                                │
                                ▼
                     [4. Mutation Detection Engine]
                     (Coordinate mapping, e.g. T790M, L858R)
                                │
                                ▼
                     [5. Distance & Phylogenetic Tree]
                     (Neighbor-Joining / UPGMA, Newick output)
                                │
                                ▼
                     [6. Visualisation & Scientific Reporting]
                     (GC skew plots, trees, mutation CSV tables)
`

## 9. Installation Instructions
`ash
git clone https://github.com/your-username/computational-bioinformatics-pipeline.git
cd computational-bioinformatics-pipeline

python -m venv venv
# On Windows:
venv\Scripts\activate
# On Linux/macOS:
source venv/bin/activate

pip install -r requirements.txt
`

## 10. Execution Instructions
### Complete Pipeline Execution
`ash
python -m src.pipeline --variants data/egfr_variants.fasta --homologs data/tp53_homologs.fasta --output results
`

### Reproducible Demonstration Script
`ash
python examples/run_bioinformatics_demo.py
`

## 11. Sample Results
### Identified Oncogenic Mutations in EGFR Kinase Domain
| Sequence ID | Alignment Position | Reference Residue | Variant Residue | Mutation Nomenclature | Clinical Significance |
| :--- | :--- | :--- | :--- | :--- | :--- |
| EGFR_T790M | 790 | Threonine (T) | Methionine (M) | p.Thr790Met | Gatekeeper resistance to 1st/2nd Gen TKIs |
| EGFR_L858R | 858 | Leucine (L) | Arginine (R) | p.Leu858Arg | Activating kinase domain mutation |

### Physicochemical Summary of EGFR Kinase Domain (Wild Type)
- **Length:** 322 amino acids
- **Molecular Weight:** 36,412.8 Da
- **Isoelectric Point (pI):** 6.42
- **GRAVY Hydrophobicity:** -0.194 (Moderately hydrophilic, globular cytoplasmic domain)
- **Instability Index:** 37.1 (Predicted stable in vitro, threshold < 40)

## 12. Screenshots & Scientific Graphs
Visualisation outputs generated in 
esults/figures/:
1. gc_skew_profile.png: Sliding-window GC skew across the nucleotide sequence highlighting local strand compositional bias.
2. phylogenetic_tree.png: Publication-grade cladogram of the TP53 gene family displaying mammalian evolutionary branch distances.

## 13. Testing Information
Automated unit tests verify IUPAC sequence validation, GC skew calculation, in-silico translation, ProtParam metrics, pairwise alignment score matrices, and mutation calling coordinates.
`ash
pytest tests/ -v
`
Test results: **9 passed, 0 failed (100% pass rate)**.

## 14. Limitations
1. **Primary Structure Assumption:** Alignments model primary residue homology without direct incorporation of crystallographic 3D tertiary contact maps.
2. **Distance-Based Phylogenetics:** Neighbor-Joining provides rapid phylogenetic reconstruction but does not model maximum likelihood evolutionary substitution rates across heterogeneous codon sites.
3. **Biological Validation:** In silico variant calling identifies nucleotide differences; clinical resistance phenotypes require verification in cell viability assays and clinical trials.

## 15. References
1. Cock, P. J., et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25(11), 1422-1423.
2. Needleman, S. B., & Wunsch, C. D. (1970). A general method applicable to the search for similarities in the amino acid sequence of two proteins. *J Mol Biol*, 48(3), 443-453.
3. Saitou, N., & Nei, M. (1987). The neighbor-joining method: a new method for reconstructing phylogenetic trees. *Mol Biol Evol*, 4(4), 406-425.

## 16. Citation Information
`ibtex
@misc{computational_bioinformatics_pipeline_2026,
  author = {Pharmaceutical Technology & Computational Research Portfolio},
  title = {Computational Bioinformatics Pipeline: Sequence Analysis, Mutation Detection, and Phylogenetic Reconstruction},
  year = {2026},
  publisher = {GitHub},
  howpublished = {\url{https://github.com/your-username/computational-bioinformatics-pipeline}}
}
`
