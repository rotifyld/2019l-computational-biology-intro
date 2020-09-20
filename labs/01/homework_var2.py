from datetime import datetime

import numpy as np

from Bio import SeqIO


def kmers(s, k):
    length = len(s)
    seq_str = s.seq._data
    return set([seq_str[i:i + k] for i in range(length - k)])


# def comp_sets(s1, s2) -> bool:

def num_to_base_2(number):
    bin_str = (bin(number)[2:])
    max_pow = len(bin_str) - 1
    pows = [int(char) * i for i, char in enumerate(reversed(bin_str))]
    bases = [int(char) * (2 ** i) for i, char in enumerate(reversed(bin_str))]
    filtered_pows = list(reversed(list(filter(lambda n: n != 0, pows))))
    if number % 2 == 1: filtered_pows.append(0)
    filtered_bases = list(reversed(list(filter(lambda n: n != 0, bases))))
    return filtered_bases, filtered_pows, max_pow


def build_full_sets(kmers_sets):
    num_sequences = len(kmers_sets)
    bases, pows, max_pow = num_to_base_2(num_sequences - 1)
    bases_cumsum = np.cumsum(bases)

    sets = np.empty((max_pow + 1, num_sequences), dtype=object)
    sets[0, :] = kmers_sets

    # go forward
    for i in range(max_pow):
        sets[i + 1, :] = sets[i, :] | np.roll(sets[i, :], -int(2 ** i))
        pass

    # go `backwards`
    for b, p, cum in zip(bases[1:], pows[1:], bases_cumsum):
        sets[max_pow, :] |= np.roll(sets[p, :], -cum)
        pass

    return sets[max_pow, :]


def compare_kmer_sets(kmers_i, kmers_but_i):
    for i, but_i in zip(kmers_i, kmers_but_i):
        if i <= but_i:
            return False
    return True


def find_smallest_k(sequences):
    sequences = list(sequences)
    shortest_sequence_length = min(map(len, sequences))

    for k in range(1, shortest_sequence_length):
        t0 = datetime.now()
        kmer_sets = [kmers(s, k) for s in sequences]
        kmer_hash_set = [set([hash(kmer) for kmer in s]) for s in kmer_sets]

        t1 = datetime.now()
        all_but_i_kmer_sets = build_full_sets(kmer_sets)

        t2 = datetime.now()
        ret = compare_kmer_sets(kmer_sets, all_but_i_kmer_sets)

        t3 = datetime.now()
        dt = t3 - t0
        print("k = {:4} kmers sets {:.2} full sets {:.2} comparing {:.2}"
              .format(k, (t1 - t0) / dt, (t2 - t1) / dt, (t3 - t2) / dt))

        if ret:
            return k


def test():
    sequences = SeqIO.parse("test_fasta.txt", "fasta")
    print(find_smallest_k(sequences))


def yeast():
    sequences = SeqIO.parse("yeast.fa", "fasta")
    print(find_smallest_k(sequences))


def test_fun():
    for k in range(3, 20):
        build_full_sets([{i} for i in range(k)])

    print('')


if __name__ == '__main__':
    # test_fun()
    # test()
    yeast()

#     k = 170, i = 175
