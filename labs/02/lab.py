from math import exp

import numpy as np
from Bio.SubsMat import SubstitutionMatrix


def substitution_matrix_to_ndarray(m: SubstitutionMatrix) -> np.ndarray:
    array = np.zeros((4, 4))
    d = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for (a, z), value in m.items():
        array[d[a], d[z]] = value
        array[d[z], d[a]] = value

    return array


def mutate(seq: str, t_step: int, q_matrix: SubstitutionMatrix, num_steps: int) -> str:
    p_delta_t = (q_matrix * t_step)
    p_delta_t_arr = substitution_matrix_to_ndarray(p_delta_t)

    # np.linalg.

    print('')


def construct_kimura_matrix(k: float) -> SubstitutionMatrix:
    data = {('A', 'A'): -2 - k, ('A', 'C'): 1, ('A', 'G'): k, ('A', 'T'): 1,
            ('C', 'C'): -2 - k, ('C', 'G'): 1, ('C', 'T'): k,
            ('G', 'G'): -2 - k, ('G', 'T'): 1,
            ('T', 'T'): -2 - k}

    return SubstitutionMatrix(data=data)


if __name__ == '__main__':
    q = construct_kimura_matrix(0.1)
    mutate('AAAAAAAAAAAAAA', 10, q, 1)
    print('')
