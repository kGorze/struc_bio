"""
we will do the needlemann_wunsch in the DP.
should use the pybenchmark to compare the performance of the code in each iteration.
make sure that the result could be in each instation. we will make it without the baseline
do it in the DP for the NW/SW. the kernel is correct for the diagonals with wavefront. kernelisation is only needed when the matrix has a lot of cells. python + numba(CPU). we should use the numba instread of the C++ compiler. this could be used when the matrice have more that 10^7 cells it will create the 10-50x than other code.
should use the Numba CPU to make it in the antydiagonals and for the middle use the prange.
numba CUDA/CuPy raw kernel - the thread is for each cell. the kernel counts it in the parellel. (NW: max{diag,up,left})
stit the antidiagonals with the wavefront. the left, up and up-left is perfect fot the DP like the whole antidiagonals k=2..(n+m) it could be parallel.

the most important idea is to how many passes should be done to find the solution. only two like (GPU,CPU) and after that (CPU).
do the Hirschberg and the Affine gap(gotoh)
what should be the traceback in the NW?

1. fill the wavefront anitidiagonal
2. traceback on the CPU and generate many paths
"""

import numpy as np
import typing
from read_fasta import read_fasta

class ResultObject:
    def __init__(self, align1: str, align2: str, score: int = 0):
        self.align1 = align1
        self.align2 = align2
        self.score = score

    def change_score(self, score: int):
        self.score = score


def print_matrix(a, n, m):
    for i in range(n):
        for j in range(m):
            print(a[i, j], end=" ")
        print()
    print()


def needleman_wunsch(seq1: str,
                     seq2: str,
                     gap: int = -2,
                     match: int = 1,
                     mismatch: int = -1) -> np.ndarray:
    seq1_len = len(seq1) + 1
    seq2_len = len(seq2) + 1

    array = np.zeros((seq2_len, seq1_len), dtype=int)

    for i in range(seq2_len):
        array[i, 0] = gap * i
    for j in range(seq1_len):
        array[0, j] = gap * j

    for i in range(1, seq2_len):
        for j in range(1, seq1_len):
            if seq2[i - 1] == seq1[j - 1]:
                score_diag = array[i - 1, j - 1] + match
            else:
                score_diag = array[i - 1, j - 1] + mismatch

            score_up = array[i - 1, j] + gap
            score_left = array[i, j - 1] + gap
            array[i, j] = max(score_diag, score_up, score_left)

    return array


def align_sequences(seq1: str,
                    seq2: str,
                    array: np.ndarray,
                    match: int,
                    mismatch: int,
                    gap: int) -> typing.List[ResultObject]:
    n = len(seq1)
    m = len(seq2)
    results: typing.List[ResultObject] = []
    final_score = array[m, n]
    dfs(n, m, "", "", results, seq1, seq2, array, match, mismatch, gap, final_score)
    return results


def dfs(i: int,
        j: int,
        align1: str,
        align2: str,
        results: typing.List[ResultObject],
        seq1: str,
        seq2: str,
        array: np.ndarray,
        match: int,
        mismatch: int,
        gap: int,
        final_score: int):
    if i == 0 and j == 0:
        res = ResultObject(align1[::-1], align2[::-1], final_score)
        results.append(res)
        return

    curr_score = array[j, i]

    if i > 0 and j > 0:
        score_diag = array[j - 1, i - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
        if curr_score == score_diag:
            dfs(i - 1, j - 1, align1 + seq1[i - 1], align2 + seq2[j - 1], results, seq1, seq2, array, match, mismatch, gap, final_score)

    if i > 0:
        score_left = array[j, i - 1] + gap
        if curr_score == score_left:
            dfs(i - 1, j, align1 + seq1[i - 1], align2 + "-", results, seq1, seq2, array, match, mismatch, gap, final_score)

    if j > 0:
        score_up = array[j - 1, i] + gap
        if curr_score == score_up:
            dfs(i, j - 1, align1 + "-", align2 + seq2[j - 1], results, seq1, seq2, array, match, mismatch, gap, final_score)



gap = -2
match = 1
mismatch = -1

sequences_1 = read_fasta("rcsb_pdb_8V8S.fasta")
sequences_2 = read_fasta("rcsb_pdb_8V91.fasta")

seq1 = sequences_1[0]
seq2 = sequences_2[0]

def calcualte_pair_file():
    pass
def calculate_from_files():
    pass
def choose_from_many_sequqnces():
    pass

def save_scores_to_file():
    with open("scores.txt", "w") as f:
        for r in best_alignments:
            f.write(f"score: {r.score}\n")
            f.write(f"{r.align1}\n")
            f.write(f"{r.align2}\n")
            f.write("\n")

def print_scores():
    best_score = max(r.score for r in all_alignments)
    best_alignments = [r for r in all_alignments if r.score == best_score]

    print(f"Number of optimal alignments: {len(all_alignments)}")
    for r in best_alignments:
        print("score:", r.score)
        print(r.align1)
        print(r.align2)
        print()

array = needleman_wunsch(seq1, seq2,
                         gap=gap,
                         match=match,
                         mismatch=mismatch)

all_alignments = align_sequences(seq1, seq2, array,
                                 match=match,
                                 mismatch=mismatch,
                                 gap=gap)

print_matrix(array, array.shape[0], array.shape[1])







if __name__ == "__main__":
    if(len(args) == 2):
        sequences_1 = read_fasta(f"{arg1}")
        sequences_2 = read_fasta(f"{arg2}")
        if(len(sequences_1) > 1):
            seq_1 = choose_from_many_sequqnces(sequences_1)
        if(len(sequences_2) > 1):
            seq_2 = choose_from_many_sequqnces(sequences_2)
        calcualte_pair_file(seq1,seq2)
    if(len(args) == 1):
        calculate_from_pair()
