from datetime import datetime
from queue import Queue

from Bio import SeqIO
from suffix_trees import STree


def bfs_find_shortest_uncommon_substring(tree: STree.STree, num_sequences):
    Q = Queue()

    curr_shortest_word = [float('inf') for _ in range(num_sequences)]
    Q.put(tree.root)

    while not Q.empty():
        node = Q.get()
        for (child, char) in node.transition_links:
            if char not in "ACGT":
                continue
            if len(child.generalized_idxs) == 1:
                idx = list(child.generalized_idxs)[0]
                word_len = node.depth + 1
                if word_len < curr_shortest_word[idx]:
                    curr_shortest_word[idx] = word_len
            else:
                Q.put(child)

    return max(curr_shortest_word)


def smallest_k(sequences):
    t0 = datetime.now()
    print(t0.time(), " experiment start")
    sequences = [s.seq._data for s in sequences]

    print(datetime.now() - t0, " building tree...")
    tree = STree.STree(sequences)
    print(datetime.now() - t0, " finished building tree, bfs...")
    ret =  bfs_find_shortest_uncommon_substring(tree, len(sequences))
    print(datetime.now() - t0, " finished all!")
    return ret


def test():
    sequences = SeqIO.parse("test_fasta.txt", "fasta")
    print(smallest_k(sequences))


def yeast():
    sequences = SeqIO.parse("yeast.fa", "fasta")
    print(smallest_k(sequences))


def tree_test():
    a = ["rpazu", "zupa"]
    tree = STree.STree(a)


if __name__ == '__main__':
    # tree_test()
    # test()
    yeast()

#
# # Generalized Suffix-Tree example.
# a = ["xxxabcxxx", "adsaabc", "ytysabcrew", "qqqabcqw", "aaabc"]
# st = STree.STree(a)
# print(st.lcs()) # "abc"
