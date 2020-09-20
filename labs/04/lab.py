from typing import List, Generator, Tuple

import graphviz
import pylab

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo import draw_ascii, draw, draw_graphviz
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor


def ex01():
    alignments = get_alignments()
    calculator = DistanceCalculator('blosum62')
    constructor = DistanceTreeConstructor()

    for a, name in alignments:
        dist_matrix = calculator.get_distance(a)
        upgma_tree = constructor.upgma(dist_matrix)
        nj_tree = constructor.nj(dist_matrix)

        print("\n\n>>> {}".format(name))
        # print(dist_matrix)
        # draw_ascii(upgma_tree)
        # draw_ascii(nj_tree)
        draw(upgma_tree)
        draw(nj_tree)


def get_alignments() -> List[Tuple[MultipleSeqAlignment, str]]:
    filenames = ['Human_PAH_paralogues.nex', 'Human_H2BFS_paralogues.nex', 'Human_PAH_orthologues.nex']

    return [(AlignIO.read(open(filename), 'nexus'), filename) for filename in filenames]


if __name__ == '__main__':
    ex01()
