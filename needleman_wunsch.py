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
2. traceback on the CPU and generate many paths\
"""

# TODO: use the pybenchmark to compare the performance of the code in each iteration of the DP. didn't have time for the rest right now. maybe I will come back to it later.

import numpy as np
import sys
import argparse
from read_fasta import read_fasta

gap = -2
match = 1
mismatch = -1


class result_object:
    def __init__(self, align1, align2, score=0):
        self.align1 = align1
        self.align2 = align2
        self.score = score

    def change_score(self, score):
        self.score = score


def print_matrix(a, n, m):
    for i in range(n):
        for j in range(m):
            print(a[i, j], "\n")
        print("\n")
    print("\n")


def needleman_wunsch(seq1, seq2, gap=-2, match=1, mismatch=-1):
    seq1_len = len(seq1) + 1
    seq2_len = len(seq2) + 1
    array = np.zeros((seq2_len, seq1_len))
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


def align_sequences(seq1, seq2, array, match, mismatch, gap, max_alignments=None):
    n = len(seq1)
    m = len(seq2)
    results = []
    final_score = array[m, n]
    if max_alignments:
        max_to_pass = max_alignments
    else:
        max_to_pass = sys.maxsize
    dfs(n, m, "", "", results, seq1, seq2, array, match, mismatch, gap, final_score, max_to_pass)
    return results


def dfs(i, j, align1, align2, results, seq1, seq2, array, match, mismatch, gap, final_score, max_alignments=10000):
    if len(results) >= max_alignments: return
    result_seq1 = align1[::-1]
    result_seq2 = align2[::-1]
    res = result_object(result_seq1, result_seq2, final_score)
    if i == 0 and j == 0:
        results.append(res)
        return
    curr_score = array[j, i]
    if i > 0 and j > 0:
        if seq1[i - 1] == seq2[j - 1]:
            score_diag = array[j - 1, i - 1] + match
        else:
            score_diag = array[j - 1, i - 1] + mismatch
        if curr_score == score_diag:
            dfs(i - 1, j - 1, align1 + seq1[i - 1], align2 + seq2[j - 1], results, seq1, seq2, array, match, mismatch,
                gap, final_score, max_alignments)
        if len(results) >= max_alignments:
            return
    if i > 0:
        score_left = array[j, i - 1] + gap
        if curr_score == score_left:
            dfs(i - 1, j, align1 + seq1[i - 1], align2 + "-", results, seq1, seq2, array, match, mismatch, gap,
                final_score, max_alignments)

    if j > 0:
        score_up = array[j - 1, i] + gap
        if curr_score == score_up:
            dfs(i, j - 1, align1 + "-", align2 + seq2[j - 1], results, seq1, seq2, array, match, mismatch, gap,
                final_score, max_alignments)


def calculate_pair_file(file):
    # this only has a file which has the two sequences in on fasta file. those are parsed and aligned.
    sequences = read_fasta(file)
    sequence_1 = sequences[0]
    sequence_2 = sequences[1]
    all_alignments = calculate_needleman_wunsch(sequence_1, sequence_2, gap, match, mismatch)
    save_scores_to_file(all_alignments)
    print_scores(all_alignments)


def calculate_from_files(file1, file2):
    sequence_1 = choose_from_many_sequences(file1)
    sequence_2 = choose_from_many_sequences(file2)
    all_alignments = calculate_needleman_wunsch(sequence_1, sequence_2, gap, match, mismatch)
    save_scores_to_file(all_alignments)
    print_scores(all_alignments)


def save_scores_to_file(all_allignments):
    with open("output.txt", "w") as f:
        for r in all_allignments:
            f.write(f"score: {r.score}\n")
            f.write(f"{r.align1}\n")
            f.write(f"{r.align2}\n\n")


def print_scores(all_alignments):
    if len(all_alignments) == 0: return
    best_score = 0
    for r in all_alignments:
        if r.score > best_score:
            best_score = r.score
    best_alignments = all_alignments[0]
    for r in all_alignments:
        if r.score == best_score:
            best_alignments = r
            break

    print(f"found: {len(all_alignments)}")
    print(f"best score: {best_score}")
    display_count = min(5, len(all_alignments))
    for idx in range(display_count):
        print(f"alignemnt {idx + 1}:")
        print(all_alignments[idx].align1)
        print(all_alignments[idx].align2, "\n")


def calculate_needleman_wunsch(seq1, seq2, gap, math, mismath, max_alignments=100):
    array = needleman_wunsch(seq1, seq2, gap=gap, match=match, mismatch=mismatch)
    all_alignments = align_sequences(seq1, seq2, array, match=match, mismatch=mismatch, gap=gap,
                                     max_alignments=max_alignments)
    return all_alignments

def choose_from_many_sequences(file):
	sequences = read_fasta(file)
	for i, seq in enumerate(sequences):
		print(f"{i}: {seq}")
	choice = int(input("choose the sequence: "))
	sequence = sequences[choice]
	return sequence


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta_file', nargs='?')
    parser.add_argument('fasta_file2', nargs='?')
    parser.add_argument('--fasta', dest='fasta_flag')
    args = parser.parse_args()
    if args.fasta_flag:
        fasta_file = args.fasta_flag
    else:
        fasta_file = args.fasta_file
    if not fasta_file: sys.exit(1)
    if args.fasta_file2:
        sequence_1 = choose_from_many_sequences(fasta_file)
        sequence_2 = choose_from_many_sequences(args.fasta_file2)
    else:
        sequences = read_fasta(fasta_file)
        if len(sequences) < 2:
            sys.exit(1)
        if len(sequences) > 2:
            sequence_1 = choose_from_many_sequences(fasta_file)
            sequence_2 = choose_from_many_sequences(fasta_file)
        else:
            sequence_1 = sequences[0]
            sequence_2 = sequences[1]
    print(f"sequence_1 len:{len(sequence_1)}")
    print(f"sequence_2 len:{len(sequence_2)}")
    all_alignments = calculate_needleman_wunsch(sequence_1, sequence_2, gap, match, mismatch)
    save_scores_to_file(all_alignments)
    print_scores(all_alignments)
