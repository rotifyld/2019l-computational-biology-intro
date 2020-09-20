from Bio import SeqIO


def kmers(s, k):
    length = len(s)
    return set([s.seq._data[i:i + k] for i in range(length - k)])

# def comp_sets(s1, s2) -> bool:
    

def exists_unique_kmers(sequences, k, prev_failed):
    num_sequences = len(sequences)
    kmer_set_per_sequence = [kmers(s, k) for s in sequences]

    unique_kmers_per_sequence = []

    r = [prev_failed]
    r.extend(list(range(num_sequences)))

    for i in r:

        all_kmers_but_i = set()
        for j in range(num_sequences):
            if i != j:
                all_kmers_but_i.update(kmer_set_per_sequence[j])

        unique_kmers = kmer_set_per_sequence[i] - all_kmers_but_i

        if len(unique_kmers) == 0:
            print("k = {}, i = {} NOT WORKING".format(k, i))
            return i
        else:
            if i % 10 == 0 or True:
                print("k = {}, i = {}    works so far".format(k, i))
            unique_kmers_per_sequence.append((i, unique_kmers))

    print(unique_kmers_per_sequence)
    return -1


def smallest_k(sequences):
    sequences = list(sequences)
    shortest_sequence_length = min(map(len, sequences))

    prev_failed_i = 0
    for k in range(shortest_sequence_length):
        ret = exists_unique_kmers(sequences, k, prev_failed_i)
        if ret == -1:
            return k
        else:
            prev_failed_i = ret
    return None


def test():
    sequences = SeqIO.parse("test_fasta.txt", "fasta")
    print(smallest_k(sequences))


def yeast():
    sequences = SeqIO.parse("yeast.fa", "fasta")
    print(smallest_k(sequences))


if __name__ == '__main__':
    test()
    # yeast()


#     k = 170, i = 175
