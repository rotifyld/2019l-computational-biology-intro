from Bio import SeqIO, AlignIO, Phylo
import Bio.Seq
from Bio.Align.Applications._Clustalw import ClustalwCommandline


def read_sequences():
    seqs_histones = list([s.seq._data for s in SeqIO.parse("histones.fa", "fasta")])


def ex01():
    seqs_histones = list(SeqIO.parse("histones.fa", "fasta"))

    for i in range(len(seqs_histones)):
        seqs_histones[i].seq = seqs_histones[i].seq.translate().rstrip("*")

    SeqIO.write(seqs_histones, "histones_trans.fa", "fasta")


def ex02():
    clustalw_cline = ClustalwCommandline("clustalw", infile="histones_trans.fa")
    clustalw_cline()
    alignment = AlignIO.read("histones_trans.aln", format="clustal")
    tree = Phylo.read("histones_trans.dnd", format="newick")

    print("")


if __name__ == '__main__':
    # ex01()
    ex02()
