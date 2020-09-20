import csv
import os
from typing import List, Iterable

from Bio.motifs.matrix import PositionSpecificScoringMatrix
from scipy.stats import binom_test
from Bio import SeqIO, motifs
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbitblastnCommandline
from Bio.Blast.Record import Blast
from Bio.motifs import Motif


def make_blast_database(in_filename: str, db_filename: str):
    """
    Make local BLAST database from given file.
    """
    NcbimakeblastdbCommandline(input_file=in_filename, parse_seqids=True, title='e_coli_genome', dbtype='nucl',
                               out=db_filename)()


def perform_local_blast_search(query_filename: str, db_filename: str, out_filename: str):
    """
    Perform local BLAST (tblastn) search and parse output to .csv file.
    """
    tmp_filename = 'data/tmp.xml'

    # perform search
    NcbitblastnCommandline(query=query_filename, db=db_filename, outfmt=5, out=tmp_filename, max_target_seqs=1)()

    # parse results
    with open(out_filename, 'w', newline='') as csv_file:
        blast_records: List[Blast] = NCBIXML.parse(open(tmp_filename))
        for record in blast_records:
            alignment = record.alignments[0]
            hsp = alignment.hsps[0]
            csv_file.write(f'{record.query[:-17]},{alignment.hit_id},{hsp.expect}\n')

    os.remove(tmp_filename)


def part01():
    """
    Delegate tasks to solve 1st part of assignment.
    """
    make_blast_database(e_coli_filename, db_filename)
    perform_local_blast_search(query_filename, db_filename, blast_out_filename)


def filter_promoters(blast_result_filename: str, promoters_filename: str,
                     groupA_promoters_filename: str, groupB_promoters_filename: str):
    """
    Filter promoters into two groups based on BLAST search result.
    """
    all_promoters = {seq.id: seq for seq in SeqIO.parse(promoters_filename, 'fasta')}

    groupA_promoters_ids = set()
    groupB_promoters_ids = set()

    with open(blast_result_filename, 'r') as f:
        reader = csv.reader(f, delimiter=',')
        for [protein_id, e_coli_id, _] in reader:
            if 'groupA' in protein_id:
                groupA_promoters_ids.add(e_coli_id)
            elif 'groupB' in protein_id:
                groupB_promoters_ids.add(e_coli_id)

    SeqIO.write((all_promoters[id] for id in groupA_promoters_ids), groupA_promoters_filename, 'fasta')
    SeqIO.write((all_promoters[id] for id in groupB_promoters_ids), groupB_promoters_filename, 'fasta')


def perform_meme_search(promoters_filename: str, out_folder: str, num_motifs: int, motif_length: int):
    """
    Prepare and execute MEME search.
    """
    command = f'meme {promoters_filename} -dna -oc {out_folder} -nostatus -nmotifs {num_motifs} -w {motif_length}'
    os.system(command)


def process_meme_output(meme_out_folder, pfm_filename):
    """
    Parse MEME search result into JASPAR .pfm files.
    """
    with open(f'{meme_out_folder}/meme.txt') as meme:
        ms = motifs.parse(meme, 'MINIMAL')
        with open(pfm_filename, 'w+') as pfm:
            pfm.write(motifs.write(ms, 'jaspar'))


def part02():
    """
    Delegate tasks to solve 2nd part of assignment.
    """

    filter_promoters(blast_out_filename, e_coli_promoters_filename,
                     e_coli_promoters_A_filename, e_coli_promoters_B_filename)

    perform_meme_search(e_coli_promoters_A_filename, meme_A_out_folder, num_motifs=10, motif_length=15)
    perform_meme_search(e_coli_promoters_B_filename, meme_B_out_folder, num_motifs=10, motif_length=15)

    process_meme_output(meme_A_out_folder, pfm_A_filename)
    process_meme_output(meme_B_out_folder, pfm_B_filename)
    pass


def parse_motifs(pfm_filename: str) -> List[Motif]:
    """
    Parse motifs from file.
    """
    with open(pfm_filename) as pfmA:
        return motifs.parse(pfmA, 'jaspar')


def get_promoter_sequences(promoters_filename) -> List[str]:
    """
    Parse sequences from file into list of strings.
    """
    with open(promoters_filename) as promoters:
        return list(record.seq._data for record in SeqIO.parse(promoters, 'fasta'))


def all_windows(s: str, window_size: int) -> Iterable[str]:
    """
    Return generator of all possible windows of given string.
    """
    return (s[i:i + window_size] for i in range(len(s) - window_size + 1))


def compute_log_odds(matrix: PositionSpecificScoringMatrix, seq: str) -> int:
    """
    Compute log odds score for given PSSM and sequence.
    """
    return sum(matrix[c][i] for i, c in enumerate(seq))


def count_hits(matrix: PositionSpecificScoringMatrix, promoters: List[str]) -> int:
    """
    Count number of hits (where log odds score > 0) for given PSSM and list of sequences.
    """
    total = 0

    for seq in promoters:
        log_odds_per_window = (compute_log_odds(matrix, window) > 0 for window in all_windows(seq, 15))
        total += sum(1 if log_odds > 0 else 0 for log_odds in log_odds_per_window)

    return total


def count_trials(seqs: List[str], window_size: int) -> int:
    """
    Count number of all possible windows of given size for given list of sequences.
    """
    return sum(len(seq) - window_size + 1 for seq in seqs)


def compute_enrichments(motifs, promoters_A, promoters_B, result_filename):
    """
    Compute enrichment scores and parse results to .csv file.
    """
    trials_A = count_trials(promoters_A, window_size=15)
    trials_B = count_trials(promoters_B, window_size=15)

    with open(result_filename, 'w', newline='') as csv_file:
        csv_file.write(
            f'motif (consensus),number of hits (A),number of hits (B),enrichment p-value (A),enrichment p-value (B)\n')

        for i, motif in enumerate(motifs):
            matrix = motif.pssm

            hits_A = count_hits(matrix, promoters_A)
            hits_B = count_hits(matrix, promoters_B)

            csv_file.write(f'{motif.consensus},{hits_A},{hits_B},'
                           f'{binom_test(hits_A, trials_A, hits_B / trials_B)},'
                           f'{binom_test(hits_B, trials_B, hits_A / trials_A)}\n')


def part03():
    """
    Delegate tasks to solve 3rd part of the assignment.
    """

    motifsA = parse_motifs(pfm_A_filename)
    motifsB = parse_motifs(pfm_B_filename)

    promoters_A = get_promoter_sequences(e_coli_promoters_A_filename)
    promoters_B = get_promoter_sequences(e_coli_promoters_B_filename)

    compute_enrichments(motifsA + motifsB, promoters_A, promoters_B, final_result_filename)


if __name__ == '__main__':
    # Dane wejściowe
    e_coli_filename = 'data/e_coli.fa'
    e_coli_promoters_filename = 'data/proms_e_coli.fa'
    query_filename = 'data/query.fa'

    # part01
    db_filename = 'data/db/db'
    blast_out_filename = 'data/blast_result.csv'

    # part02
    e_coli_promoters_A_filename = 'data/proms_e_coli_A.fa'
    e_coli_promoters_B_filename = 'data/proms_e_coli_B.fa'
    meme_A_out_folder = 'data/memeA'
    meme_B_out_folder = 'data/memeB'
    pfm_A_filename = 'data/A.pfm'
    pfm_B_filename = 'data/B.pfm'

    # part03
    final_result_filename = 'data/final_result.csv'

    part01()
    part02()
    part03()
