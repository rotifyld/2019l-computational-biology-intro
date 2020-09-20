from datetime import datetime
import logging

import numpy as np
from Bio import SeqIO


def kmers(s, k):
    length = len(s)
    return set([s[i:i + k] for i in range(length - k + 1)])


def split_by_uniqueness(sequences, unique_indices, not_unique_indices, k):
    kmers_per_sequence = [kmers(s, k) for s in sequences]

    kmers_checked = set()
    for i in unique_indices:
        kmers_checked.update(kmers_per_sequence[i])

    not_unique = set()
    unique = set()

    for i in not_unique_indices:
        i_kmers = kmers_per_sequence[i]

        all_but_i_kmers = set(kmers_checked)

        for j in not_unique_indices:
            if j > i:
                all_but_i_kmers.update(kmers_per_sequence[j])

        if i_kmers.issubset(all_but_i_kmers):
            not_unique.add(i)
        else:
            unique.add(i)

        kmers_checked.update(i_kmers)

    print("   unique = {}\n   not_unique = {}".format(unique, not_unique))

    return unique, not_unique


def binsearch_k(sequences):
    num_sequences = len(sequences)
    shortest_sequence_length = min(map(len, sequences))

    low = 1
    high = shortest_sequence_length  # exclusive (upper bound)

    # list of indices such that corresponding sequence *does not have* unique k-mer of length lower_bound - 1
    not_unique_indices = set(range(num_sequences))
    # list of indices such that corresponding sequence *has* unique k-mer of length lower_bound - 1
    unique_indices = set()

    # benchmarking
    t0 = datetime.now()
    t_last = t0

    # bin-search
    while low < high:
        mid = low
        print("new mid = {}".format(mid))

        unique, not_unique = split_by_uniqueness(sequences, unique_indices, not_unique_indices, mid)
        if len(not_unique) == 0:
            high = low
        else:
            low = mid + 1
            unique_indices.update(unique)
            not_unique_indices.difference_update(unique)

        t_curr = datetime.now()
        str = "mid = {}, len_not_unique = {}, time_single = {}, time_total = {}" \
            .format(mid, len(not_unique_indices), t_curr - t_last, t_curr - t0)
        print(str)
        logging.debug(str)
        t_last = t_curr

    assert low == high
    return low


def test():
    sequences = SeqIO.parse("test_fasta.txt", "fasta")
    sequences = [s.seq._data for s in sequences]
    print(binsearch_k(sequences))


def yeast():
    sequences = SeqIO.parse("yeast.fa", "fasta")
    sequences = [s.seq._data for s in sequences]
    print(binsearch_k(sequences))


if __name__ == '__main__':
    logging.basicConfig(filename='{}_homework_var3.log'.format(datetime.now()), level=logging.DEBUG)

    # test()
    yeast()
