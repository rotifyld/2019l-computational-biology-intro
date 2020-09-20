from typing import Tuple, List

from Bio import SeqIO, Seq
from Bio.pairwise2 import align, format_alignment
from Bio.SubsMat import MatrixInfo

# sequence of nucleotides (over alphabet Σ = {A, C, G, T}).
SequenceN = str

# sequence of amino acids.
SequenceAA = str

# local alignment, as returned by Bio.pairwise2.alignXX function
LocalAlignment = Tuple[str, str, float, int, int]


class SequenceTranslation:
    """
    Sequence represented as both nucleotides and amino acids.
    """

    def __init__(self, seq_n: SequenceN, seq_aa: SequenceAA = None):
        self.seq_n = seq_n
        self.seq_aa = Seq.translate(seq_n) if seq_aa is None else seq_aa

    def split(self) -> List['SequenceTranslation']:
        """
        List of sequences split by "*" char in amino acidic representation.
        In special case of original sequence not having "*" char, returns a singleton of itself.
        """
        seq_aa = self.seq_aa
        seq_n = self.seq_n
        ret = []
        pos = None

        while pos != -1:
            pos = seq_aa.find("*")

            if pos == -1:
                if len(seq_aa) > 0:
                    new_seq = SequenceTranslation(seq_n, seq_aa)
                    ret.append(new_seq)
            else:
                new_seq_aa = seq_aa[:pos]
                new_seq_n = seq_n[:3 * pos]
                if len(new_seq_aa) > 0:
                    new_seq = SequenceTranslation(new_seq_n, new_seq_aa)
                    ret.append(new_seq)
                seq_aa = seq_aa[pos + 1:]
                seq_n = seq_n[3 * (pos + 1):]

        return ret


def possible_translations(seq: SequenceN) -> List[SequenceTranslation]:
    """
    :return: list of all (3) possible translations, by discarding first 0, 1, or 2 nucleotides.
    """

    def truncate_tail(s: SequenceN) -> SequenceN:
        """ Helper function returning truncated sequence so that len(seq) mod 3 == 0. """
        return s if len(s) % 3 == 0 else s[:-(len(s) % 3)]

    return [SequenceTranslation(truncate_tail(seq[i:])) for i in range(3)]


class LocalAlignmentWrapper:
    """
    Wrapper class for alignment: computing, comparing, and converting to nucleotides.
    """

    def __init__(self, seq1: SequenceTranslation, seq2: SequenceTranslation,
                 match_dict: dict, open_gap_penalty: float, extend_gap_penalty: float):
        self.seq1 = seq1
        self.seq2 = seq2

        alignments = align.localds(seq1.seq_aa, seq2.seq_aa, match_dict, open_gap_penalty, extend_gap_penalty,
                                   one_alignment_only=True)

        if len(alignments) == 0:
            self.value = float('-inf')
        else:
            _, _, self.value, _, _ = alignments[0]
            self.alignment_aa = alignments[0]

    def __lt__(self, other):
        if not isinstance(other, LocalAlignmentWrapper):
            return False
        else:
            return self.value < other.value

    def alignment_n(self) -> LocalAlignment:
        """
        Returns the alignment but with amino acids converted back to nucleotides.
        """

        def aligned_aa_to_n(seq_n: SequenceN, seq_aa_aligned: SequenceAA) -> SequenceN:
            """ Helper function returning seq_n but extended with gaps in places where seq_aa_aligned has gaps. """
            tmp = []
            non_gap = 0
            for i, c in enumerate(seq_aa_aligned):
                if c == '-':
                    tmp.append('---')
                else:
                    tmp.append(seq_n[3 * non_gap:3 * non_gap + 3])
                    non_gap += 1

            return "".join(tmp)

        seq1_aligned_aa, seq2_aligned_aa, value, lower_bound, upper_bound = self.alignment_aa

        return (aligned_aa_to_n(self.seq1.seq_n, seq1_aligned_aa), aligned_aa_to_n(self.seq2.seq_n, seq2_aligned_aa),
                value, 3 * lower_bound, 3 * upper_bound)


def best_local_alignment(seq1_n: SequenceN, seq2_n: SequenceN,
                         open_gap_penalty: float, extend_gap_penalty: float) -> LocalAlignment:
    """
    Finds best local alignment as required in task specification.
    """

    # translate sequences into amino acids (in all 3 possible ways) and split them by stop codon
    seq1_all_fragments = [subseq_t for seq_t in possible_translations(seq1_n) for subseq_t in seq_t.split()]
    seq2_all_fragments = [subseq_t for seq_t in possible_translations(seq2_n) for subseq_t in seq_t.split()]

    alignments = []
    for s1 in seq1_all_fragments:
        for s2 in seq2_all_fragments:
            alignment = LocalAlignmentWrapper(s1, s2,
                                              match_dict=MatrixInfo.blosum60,
                                              open_gap_penalty=open_gap_penalty,
                                              extend_gap_penalty=extend_gap_penalty)
            alignments.append(alignment)

    best_alignment = max(alignments)
    return best_alignment.alignment_n()


def histones():
    seqs_histones = list([s.seq._data for s in SeqIO.parse("histones.fa", "fasta")])

    for i, s1 in enumerate(seqs_histones):
        for j, s2 in enumerate(seqs_histones):
            if i < j:
                alignment = best_local_alignment(s1, s2, open_gap_penalty=-1, extend_gap_penalty=-0.5)
                print(format_alignment(*alignment))


if __name__ == '__main__':
    histones()
