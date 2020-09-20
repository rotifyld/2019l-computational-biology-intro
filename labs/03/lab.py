from typing import Tuple

from Bio import Seq, SeqIO, SubsMat, pairwise2
from Bio.SubsMat import MatrixInfo


def compute_pair_ex3(s1: str, s2: str) -> Tuple[float, float]:
    return (pairwise2.align.localmd(s1, s2, 1, -1, -1, -0.5, -1, -0.5, score_only=True),
            pairwise2.align.globalmd(s1, s2, 1, -1, -1, -0.5, -1, -0.5, score_only=True))


blosum60 = MatrixInfo.blosum60


def compute_pair_ex4(s1: str, s2: str) -> Tuple[float, float]:
    return (pairwise2.align.localdd(s1, s2, blosum60, -1, -0.5, -1, -0.5, score_only=True),
            pairwise2.align.globaldd(s1, s2, blosum60, -1, -0.5, -1, -0.5, score_only=True))


def compute_all_with_all(sequences, sequences2=None, function=None):
    assert function is not None

    total_local, total_global = 0, 0
    num = 0
    for i, s1 in enumerate(sequences):
        for j, s2 in enumerate(sequences if sequences2 is None else sequences2):
            if sequences2 is not None or i < j:
                num += 1
                loc, glob = function(s1, s2)
                total_local += loc
                total_global += glob
                # print("  {:2} -- {:2}: local = {:4} global = {:4}".format(i, j, int(loc), int(glob)))

    mean_local = total_local / num
    mean_global = total_global / num
    print("mean: local = {:4} global = {:4}".format(int(mean_local), int(mean_global)))


def read_sequences():
    seqs_histones = list([s.seq._data for s in SeqIO.parse("histones.fa", "fasta")])
    seqs_bzips = list([s.seq._data for s in SeqIO.parse("bzips.fa", "fasta")])
    return seqs_histones, seqs_bzips


def ex3():
    seqs_histones, seqs_bzips = read_sequences()

    print("histones:")
    compute_all_with_all(seqs_histones, function=compute_pair_ex3)
    print("bzips:")
    compute_all_with_all(seqs_bzips, function=compute_pair_ex3)
    print("bzips x histones")
    compute_all_with_all(seqs_histones, seqs_bzips, function=compute_pair_ex3)


def ex4():
    seqs_histones, seqs_bzips = read_sequences()
    seqs_histones = [Seq.translate(s, to_stop=True) for s in seqs_histones]
    seqs_bzips = [Seq.translate(s, to_stop=True) for s in seqs_bzips]

    print("histones:")
    compute_all_with_all(seqs_histones, function=compute_pair_ex4)
    print("bzips:")
    compute_all_with_all(seqs_bzips, function=compute_pair_ex4)
    print("bzips x histones")
    compute_all_with_all(seqs_histones, seqs_bzips, function=compute_pair_ex4)


if __name__ == '__main__':
    ex3()
    # ex4()
