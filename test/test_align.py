# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1,header = read_fasta("./data/test_seq1.fa")
    seq2,header = read_fasta("./data/test_seq2.fa")

    # Matrix dimensions are correct
    # M[0][0] = 0
    # Final score in matrices matches alignment_score

    nw = NeedlemanWunsch(
        sub_matrix_file="./substitution_matrices/BLOSUM62.mat", 
        gap_open=-10,
        gap_extend=-1
    )

    score, seqA, seqB = nw.align(seq1, seq2)

    # Check matrix dimensions
    assert nw._align_matrix.shape == (len(seq1)+1, len(seq2)+1)
    assert nw._gapA_matrix.shape == (len(seq1)+1, len(seq2)+1)
    assert nw._gapB_matrix.shape == (len(seq1)+1, len(seq2)+1)

    # Check M[0][0] = 0
    assert nw._align_matrix[0][0] == 0

    # Check alignment score matches best score at bottom-right
    final_scores = [
        nw._align_matrix[len(seq1)][len(seq2)],
        nw._gapA_matrix[len(seq1)][len(seq2)],
        nw._gapB_matrix[len(seq1)][len(seq2)]
    ]
    assert score == max(final_scores)   

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3,header = read_fasta("./data/test_seq3.fa")
    seq4,header = read_fasta("./data/test_seq4.fa")
    
    # Create NW object
    # Run alignment
    # Assert the score is 18
    # Assert the alignment strings are correct

    nw = NeedlemanWunsch(
        sub_matrix_file="./substitution_matrices/BLOSUM62.mat", 
        gap_open=-10,
        gap_extend=-1
    )

    score, seqA, seqB = nw.align(seq3, seq4)
    assert score == 18
    assert seqA == "MAVHQLIRRP"
    assert seqB == "M---QLIRHP"




