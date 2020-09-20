from Bio import AlignIO, Phylo, SeqIO, Entrez
from Bio.Align.Applications import ClustalwCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from matplotlib import pyplot as plt

import os

DATA_DIR = "data"


# Download complete sequence (and optionally spike genome) of a single virus from entrez database with specified terms
def entrez_download(genome_search_term, spike_search_term, filename):
    print("Downloading {}...".format(filename))

    # complete genome
    search_handle = Entrez.esearch(db="nuccore", term=genome_search_term)
    search = Entrez.read(search_handle)
    search_handle.close()

    record_handle = Entrez.efetch(db="nucleotide", id=search["IdList"][0], rettype="fasta")
    record = SeqIO.read(record_handle, "fasta")
    record_handle.close()

    SeqIO.write(record, "{}/{}.fa".format(DATA_DIR, filename), "fasta")

    # spike protein
    if spike_search_term is not None:
        search_handle = Entrez.esearch(db="protein", term=spike_search_term)
        search = Entrez.read(search_handle)
        search_handle.close()

        spike_handle = Entrez.efetch(db="protein", id=search["IdList"][0], rettype="fasta")
        record = SeqIO.read(spike_handle, "fasta")
        spike_handle.close()

        SeqIO.write(record, "{}/{}_spike.fa".format(DATA_DIR, filename), "fasta")


# Download all viruses as specified in exercise
def download_viruses():
    entrez_download("Severe acute respiratory syndrome coronavirus 2 complete genome",
                    "Severe acute respiratory syndrome coronavirus 2 spike protein",
                    filename="SARS-CoV-2")

    entrez_download("SARS coronavirus ZJ0301 complete genome",
                    "SARS coronavirus ZJ0301 spike protein",
                    filename="SARS-CoV")

    entrez_download("Bat SARS-like coronavirus WIV1 complete genome",
                    "Bat SARS-like coronavirus WIV1 spike protein",
                    filename="Bat_SL-CoV-WIV1")

    entrez_download("MERS coronavirus complete genome",
                    "Middle East respiratory syndrome-related coronavirus isolate spike protein",
                    filename="MERS-CoV")

    entrez_download("Hepatovirus A complete genome", None, filename="Hepatitis_A")

    # Download influenza segments
    for i in range(1, 8):
        entrez_download("Influenza A A+chicken+Jiangsu+1001+2013+H5N2 segment+{}".format(i), None,
                        filename="Influenza_A_{}".format(i))

    # Concatenate Influenza segments
    influenza_concatenated = SeqIO.SeqRecord("", description="Influenza A genome concatenated segments")
    for i in range(1, 8):
        influenza_fragment = SeqIO.parse("{}/Influenza_A_{}.fa".format(DATA_DIR, i), "fasta")
        influenza_concatenated += next(influenza_fragment)
        os.remove("{}/Influenza_A_{}.fa".format(DATA_DIR, i))
    SeqIO.write(influenza_concatenated, "{}/Influenza_A.fa".format(DATA_DIR), "fasta")


def preprocess_data():
    genomes = []
    spikes = []
    for _, _, files in os.walk(DATA_DIR):
        for file in files:
            [name, extention] = file.split(".")
            if extention == 'fa':
                for seq in SeqIO.parse("{}/{}".format(DATA_DIR, file), "fasta"):
                    if name.endswith("_spike"):
                        name = name[:-6]
                        seq.id = name
                        spikes.append(seq)
                    else:
                        seq.id = name
                        genomes.append(seq)

    SeqIO.write(genomes, "complete_sequences.fa", "fasta")
    SeqIO.write(spikes, "spike_proteins.fa", "fasta")


def get_data():
    # Setup search engine
    Entrez.email = "db394094@students.mimuw.edu.pl"
    Entrez.tool = "testing"

    # Create data directory if don't exist
    try:
        os.mkdir(DATA_DIR)
    except Exception as e:
        pass

    download_viruses()
    preprocess_data()


def build_trees(filename, tree_name):
    # Compute alignment with ClustalW algorithm
    clustalw_cline = ClustalwCommandline("clustalw", infile="{}.fa".format(filename))
    clustalw_cline()
    alignment = AlignIO.read("{}.aln".format(filename), format="clustal")

    # Create distance matrix
    calculator = DistanceCalculator('blosum62')
    dist_matrix = calculator.get_distance(alignment)

    # Build phylogenetic trees using upgma and nj methods
    constructor = DistanceTreeConstructor()
    upgma_tree = constructor.upgma(dist_matrix)
    nj_tree = constructor.nj(dist_matrix)

    # Draw the trees
    label_func = lambda clade: "" if clade.name.startswith("Inner") else clade

    Phylo.draw(upgma_tree, label_func=label_func, do_show=False)
    plt.title("{} × upgma".format(tree_name))
    plt.show()

    Phylo.draw(nj_tree, label_func=label_func, do_show=False)
    plt.title("{} × nj".format(tree_name))
    plt.show()


if __name__ == '__main__':
    get_data()
    build_trees("spike_proteins", "spike proteins")
    build_trees("complete_sequences", "complete genome")
