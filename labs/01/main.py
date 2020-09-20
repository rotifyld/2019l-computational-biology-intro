from collections import deque
from typing import Iterable, Generator

from Bio import Seq, SeqRecord, SeqIO


def load_sequences():
    return SeqIO.parse("test_fasta.txt", "fasta")


def ex01():
    sequences = load_sequences()
    rc_sequences = [s.reverse_complement(id=True, description=True) for s in sequences]
    SeqIO.write(rc_sequences, "test_reverse_complementary.fasta", "fasta")


def kmers(s, k):
    length = len(s)
    return set([s._data[i:i + k] for i in range(length - k)])


def kmers_comp(s, k):
    kmers(s, k).union(kmers(s.reverse_comlementary(), k))


def ex02():
    sequences = SeqIO.parse("test_fasta.txt", "fasta")
    for s in sequences:
        for k in range(1, len(s)):
            print("seq = {}, k = {}, kmers = {}".format(s.seq, k, kmers(s.seq, k)))


def create_euler_graph(kmers):
    V = set()
    E = {}
    E_ = {}
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        V.add(prefix)
        V.add(suffix)

        prefix_out = E.get(prefix, list())
        prefix_out.append(suffix)
        E[prefix] = prefix_out

        suffix_in = E_.get(suffix, list())
        suffix_in.append(prefix)
        E_[suffix] = suffix_in

    return V, E, E_


def euler(kmers):
    kmers = set(kmers)
    V, E, E_ = create_euler_graph(kmers)

    potential_starting_points = list(filter(lambda v: len(E.get(v, [])) - len(E_.get(v, [])) == 1, V))
    print(potential_starting_points)

    for v0 in potential_starting_points:
        Q = deque()
        Q.append((v0, {v0}))

        (v, visited) = Q.pop()
        while visited != V:
            next = set(E.get(v, []))
            possible_next = next.difference(visited)
            if len(possible_next) == 0:
                return -1
            else:
                u = possible_next


def ex03():
    sequences = SeqIO.parse("test_fasta.txt", "fasta")
    for s in sequences:
        for k in range(3, len(s)):
            print("seq = {}, k = {}, euler_paths = {}".format(s.seq, k, euler(kmers(s.seq, k))))


if __name__ == '__main__':
    ex03()

# zad dom:
# optymalne rozwiązanie: tablice sufixowe O(n)
# zazwwyczaj o(n^2) i wykorzystuje hashtable lub słownik
