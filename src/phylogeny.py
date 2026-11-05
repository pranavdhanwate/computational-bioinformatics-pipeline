from io import StringIO
from typing import List, Tuple
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio.SeqRecord import SeqRecord

from .alignment import SequenceAligner


class PhylogeneticTreeBuilder:
    """Constructs phylogenetic trees using pairwise sequence distance matrices."""

    def __init__(self, method: str = "nj"):
        self.method = method.lower()
        self.aligner = SequenceAligner(mode="global", is_protein=False)
        self.constructor = DistanceTreeConstructor()

    def build_tree_from_records(self, records: List[SeqRecord]) -> Tuple[Phylo.BaseTree.Tree, str]:
        names = [rec.id.split("|")[0] for rec in records]
        n = len(names)
        matrix = []

        for i in range(n):
            row = []
            for j in range(i + 1):
                if i == j:
                    row.append(0.0)
                else:
                    aln = self.aligner.align(str(records[i].seq), str(records[j].seq))
                    distance = round(max(0.0, 1.0 - (aln["identity_percent"] / 100.0)), 4)
                    row.append(distance)
            matrix.append(row)

        dm = DistanceMatrix(names, matrix)

        if self.method == "upgma":
            tree = self.constructor.upgma(dm)
        else:
            tree = self.constructor.nj(dm)

        # Generate Newick string
        newick_io = StringIO()
        Phylo.write(tree, newick_io, "newick")
        newick_str = newick_io.getvalue().strip()

        return tree, newick_str
