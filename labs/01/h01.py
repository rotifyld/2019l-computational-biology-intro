from datetime import datetime
from typing import List, Set, Tuple

from Bio import SeqIO


def kmers(s: str, k: int) -> Set[str]:
    """
    :return: Set of all k-mers of the sequence s
    """
    length = len(s)
    return set([s[i:i + k] for i in range(length - k + 1)])


def split_by_uniqueness(sequences: List[str],
                        indices_unique: Set[int],
                        indices_not_unique: Set[int],
                        k: int) \
        -> Tuple[Set[int], Set[int]]:
    """
    Assuming `sequences` - a list of sequences that are divided into:
        - sequences having unique m-mers (set of indices of such sequences: `indices_unique`),
        - sequences not having unique m-mers (set of indices of such sequences: `indices_not_unique`),
    such that there exists integer m < `k`.

    Splits `indices_not_unique` into two subsets:
        - sequences having unique `k`-mers,
        - sequences not having unique `k`-mers.
    """
    kmers_per_sequence = [kmers(s, k) for s in sequences]

    # For performance reasons, to compute set of k-mers of all but one sequence.
    kmers_checked = set()
    for i in indices_unique:
        kmers_checked.update(kmers_per_sequence[i])

    # Returned values of split set.
    # Their intersection is always empty.
    # At the return their sum is equal to `indices_not_unique`.
    unique = set()
    not_unique = set()

    for i in indices_not_unique:

        # set of k-mers of i-th sequence.
        i_kmers = kmers_per_sequence[i]

        # construct set of k-mers of all sequences but the i-th.
        all_but_i_kmers = set(kmers_checked)
        for j in indices_not_unique:
            if j > i:
                all_but_i_kmers.update(kmers_per_sequence[j])

        # compare above two sets and update
        if i_kmers.issubset(all_but_i_kmers):
            not_unique.add(i)
        else:
            unique.add(i)

        kmers_checked.update(i_kmers)

    return unique, not_unique


def search_k(sequences: List[str]) -> int:
    """
    Performs binary search in order to find smallest k satisfying the assignment condition.
    """
    num_sequences = len(sequences)
    shortest_sequence_length = min(map(len, sequences))

    low = 1
    high = shortest_sequence_length  # exclusive (upper bound)

    # For performance reasons.
    # For meaning, see: loop invariant.
    not_unique_indices = set(range(num_sequences))
    unique_indices = set()

    # Binary search loop
    while low < high:
        # while loop invariants:
        #   There does not exist a solution for problem of finding a set of probes of length `low - 1`
        #   There exists a solution for problem of finding a set of probes of length `high`
        #
        #   `not_unique_indices` - set of indices such that sequence does not have unique k-mer of length `low - 1`
        #   `unique_indices` - set of indices such that sequence has unique k-mer of length `low - 1`

        mid = (low + high) // 2
        unique, not_unique = split_by_uniqueness(sequences, unique_indices, not_unique_indices, mid)

        if len(not_unique) == 0:
            high = mid
        else:
            low = mid + 1
            unique_indices.update(unique)
            not_unique_indices.difference_update(unique)

    return low  # at this point, `low == high`


def yeast():
    sequences = SeqIO.parse("test_fasta.txt", "fasta")
    sequences = [s.seq._data for s in sequences]  # convert List[SeqRecords] to List[srt]
    t0 = datetime.now()

    k = search_k(sequences)
    print("The solution: k = {}\nFinding it took: {}".format(k, datetime.now() - t0))


if __name__ == '__main__':
    yeast()
