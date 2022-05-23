#!/usr/bin/env python3

"""
Author: Noël Huang

Description: this is a script to align two amino acid sequences, based on the Needleman-Wunsch algorithm. It makes use
of dynamic programming, and utilizes the BLOSUM62 scoring matrix. The output is a provably optimal alignment of the
two given sequences, for a given end gap penalty and regular gap penalty. Sequences must be given inside this python
script (no parsing function is available (yet)).

Usage: py alignsequence.py
"""
# import statements here
import numpy as np

# functions between here and __main__
blosum = """
# http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
   A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
   R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
   N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
   D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
   C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
   Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
   E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
   G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
   H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
   I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
   L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
   K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
   M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
   F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
   P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
   S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
   T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
   W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
   Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
   V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
   B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
   Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
   X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
   * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 
"""


def blosum62():
    """Return order and similarity scores from BLOSUM62 matrix

    order: dict of {res: idx_in_matrix}
    blosum_matrix: list of lists with similarity scores
    """
    order = {}
    blosum_matrix = []
    for line in blosum.split('\n'):
        if line.startswith('#'):
            continue
        if not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) == 24:
            for idx, sym in enumerate(parts):
                order[sym] = idx
        else:
            # list around the map construction for python3 compatibility
            blosum_matrix.append(list(map(int, parts[1:])))
    return order, blosum_matrix


BLOSUM62_ORDER, BLOSUM62_MATRIX = blosum62()


def score(res1, res2):
    """Return similarity score from BLOSUM62 matrix for two residues
    
    res1: string, amino acid
    res2: string, amino acid
    """
    lookup1 = BLOSUM62_ORDER[res1]
    lookup2 = BLOSUM62_ORDER[res2]
    return BLOSUM62_MATRIX[lookup1][lookup2]


# write your own functions below here


def elongate_sequences(sequence):
    sequence = " " + sequence
    #print(sequence)
    return sequence


def create_matrix(seq1, seq2):
    matrix = np.array([[None for j in range(len(seq2))] for i in range(len(seq1))])
    #print(f"This is an empty matrix: \n {matrix}")
    return matrix


def traceback(filled_matrix, seq1, seq2, gap_penalty, end_gap_penalty, i, j, alignment):
    """
    Starts at starting position i, j and traces back the optimal path in the dynamic programming scoring matrix.

    Args:
        filled_matrix: numpy array: The dynamic programming matrix filled with optimal scores for each cell
        seq1: String: First amino acid sequence, ONLY capital letters
        seq2: String: Second amino acid sequence, ONLY capital letters
        gap_penalty: int: The penalty that should be given for gaps in the middle of the alignment.
                                Should be negative integer.
        end_gap_penalty: int: The penalty that should be given for gaps at the beginning and end of the alignment.
                                Should be negative integer.
        i: int: starting position in the filled_matrix (row index).
        j: int: starting position in the filled_matrix (bottom index).
        alignment: list of two strings, tail of the alignment.

    Returns:
        alignment: list of two strings, each string is a sequence of amino acids. Gaps are indicated by a dash.
        The index of one string is aligned with the index of the other string.
    """
    while i > 0 and j > 0:
        diagonal_score = score(seq1[i], seq2[j])  # Score associated with diagonal movement

        if j == len(filled_matrix[i]) - 1:  # "If we are in the last column":
            vertical_score = end_gap_penalty  # "Apply the end_gap penalty"
        else:  # "Else, apply normal gap penalty"
            vertical_score = gap_penalty

        if i == len(filled_matrix) - 1:  # "If we are in the last row":
            horizontal_score = end_gap_penalty  # "Apply the end_gap penalty"
        else:  # "Else, apply normal gap penalty"
            horizontal_score = gap_penalty

        if filled_matrix[i][j] - horizontal_score == filled_matrix[i][j - 1]:
            # A horizontal movement was performed, thus there is a gap in top alignment string
            alignment[0] = "-" + alignment[0]
            alignment[1] = seq2[j] + alignment[1]
            j = j-1
        elif filled_matrix[i][j] - diagonal_score == filled_matrix[i - 1][j - 1]:
            # A diagonal movement was performed, hence index i was aligned with index j
            alignment[0] = seq1[i] + alignment[0]
            alignment[1] = seq2[j] + alignment[1]
            i = i-1
            j = j-1
        elif filled_matrix[i][j] - vertical_score == filled_matrix[i - 1][j]:
            # A vertical movement was performed, thus there is a gap in bottom alignment string
            alignment[0] = seq1[i] + alignment[0]
            alignment[1] = "-" + alignment[1]
            i = i-1

    while i > 0:  # If the optimal path enters the leftmost column before 0,0 , we must still finish the alignment
        alignment[0] = seq1[i] + alignment[0]
        i = i-1
    while j > 0:  # If the optimal path enters the top row before 0,0 , we must still finish the alignment
        alignment[1] = seq2[j] + alignment[1]
        j = j-1

    # If the two strings in the alignment are of different length, gaps are added in front of the shorter string to make
    # it of equal length to the other string.
    if len(alignment[0]) > len(alignment[1]):
        length_difference = len(alignment[0]) - len(alignment[1])
        alignment[1] = length_difference * "-" + alignment[1]
    elif len(alignment[1]) > len(alignment[0]):
        length_difference = len(alignment[1]) - len(alignment[0])
        alignment[0] = length_difference * "-" + alignment[0]

    return alignment


def fill_matrix(seq1, seq2, matrix, end_gap_penalty, gap_penalty):
    """
    Aligns two sequences optimally, prints alignment, and returns the dynamic programming scoring matrix.

    Args:
        seq1: String: First amino acid sequence, ONLY capital letters
        seq2: String: Second amino acid sequence, ONLY capital letters
        matrix: Numpy array: A matrix of size len(seq1) by len(seq2) which contains None in all cells.
        end_gap_penalty: int: The penalty that should be given for gaps at the beginning and end of the alignment.
                                Should be negative integer.
        gap_penalty: int: The penalty that should be given for gaps in the middle of the alignment.
                                Should be negative integer.

    Returns:
        filled_matrix: The dynamic programming matrix filled with optimal scores for each cell.

    This function also prints several noteworthy data:
        The actual alignment, the percent identity, and the alignment score.
    Additionally, mentions of i and j in this function refer the the tow index (i) and column index (j)

    """
    # traceback_matrix = np.array([[0 for j in range(len(seq2))] for i in range(len(seq1))])
    # print(f"This is a traceback matrix: \n {traceback_matrix}")

    matrix[0][0] = 0  # Starting point (top-left matrix) has 0 score
    for j in range(len(matrix[0])):
        matrix[0][j] = 0
    for i in range(1, len(matrix)):
        matrix[i][0] = 0

    for j in range(1, len(matrix[0])):  # All values in top row can be calculated from end_gap_penalty
        matrix[0][j] = matrix[0][j - 1] + end_gap_penalty

    for i in range(1, len(matrix)):  # All values in leftmost column can be calculated from end_gap_penalty
        matrix[i][0] = matrix[i - 1][0] + end_gap_penalty

    for i in range(1, len(matrix)):  # Find the highest possible score for every cell in the matrix.
        for j in range(1, len(matrix[i])):

            eq1 = matrix[i - 1][j - 1] + score(seq1[i], seq2[j])  # Score for diagonal movement in matrix
            # (corresponds with aligning two amino acids)

            if j == len(matrix[i]) - 1:  # "If we are in the last column":
                eq2 = matrix[i - 1][j] + end_gap_penalty  # "Apply the end_gap penalty"
            else:  # "Else, apply normal gap penalty"
                eq2 = matrix[i - 1][j] + gap_penalty

            if i == len(matrix) - 1:  # "If we are in the last row":
                eq3 = matrix[i][j - 1] + end_gap_penalty  # "Apply the end_gap penalty"
            else:  # "Else, apply normal gap penalty"
                eq3 = matrix[i][j - 1] + gap_penalty

            possible_scores = [eq1, eq2, eq3]
            matrix[i][j] = max(possible_scores)  # Cell i j is updated with maximum possible score
    filled_matrix = matrix
    print(f"Based on: end gap penalty = {end_gap_penalty}, gap penalty = {gap_penalty}")
    print(f"This is the matrix with scores: \n {filled_matrix}")

    # Introduce "alignment", a list of two strings which will form the alignment
    # (where index of first string is aligned with index of second string)
    alignment = ["", ""]

    if end_gap_penalty < 0:  # If there is an end gap penalty, we know that optimal path ends in bottom right
        i = len(filled_matrix) - 1  # i = Row index
        j = len(filled_matrix[i]) - 1  # j = Column index
    else:  # If there is no end gap penalty, we must find the highest scoring cell at the edges and build from there.
        max_score_last_row = max(filled_matrix[-1])
        last_column = [filled_matrix[i][-1] for i in range(len(filled_matrix))]
        max_score_last_column = max(last_column)
        if max_score_last_row >= max_score_last_column:  # If the max score is in the last row, find coordinates:
            i = len(filled_matrix) - 1
            j = list(filled_matrix[-1]).index(max_score_last_row)

            # Here we make end gaps if there should be any
            number_of_gaps = len(filled_matrix[-1]) - 1 - j
            # "len(filled_matrix[-1]) - 1" denotes the highest index in the bottom row

            alignment[0] = number_of_gaps * "-" + alignment[0]

            for index in range(j+1, j + number_of_gaps+1): # Start adding amino acids to the tail of the alignment
                alignment[1] = alignment[1] + seq2[index]
        else:  # If the max score is in the last column, find coordinates:
            j = len(filled_matrix[0]) - 1
            i = last_column.index(max_score_last_column)
            number_of_gaps = len(filled_matrix) - 1 - i
            alignment[1] = number_of_gaps * "-" + alignment[1]  # Make gaps at the end if necessary

            for index in range(i+1, i + number_of_gaps+1):  # Start adding amino acids to the tail of the alignment
                alignment[0] = alignment[0] + seq1[index]

    traceback(filled_matrix, seq1, seq2, gap_penalty, end_gap_penalty, i, j, alignment)  # Start traceback function

    line_length = 60  # This block just makes prettier print output
    lines_0 = [alignment[0][i:i + line_length] + '\n' for i in range(0, len(alignment[0]), line_length)]
    lines_1 = [alignment[1][i:i + line_length] + '\n' for i in range(0, len(alignment[1]), line_length)]

    print("This is the alignment:")
    for i in range(0, len(lines_0)):
        print(f"{lines_0[i]}{lines_1[i]}")

    number_identical = 0
    for i in range(0, len(alignment[0])):
        if alignment[0][i] == alignment[1][i]:
            number_identical += 1
    percent_identity = number_identical/len(alignment[0])*100
    print(f"This is the percent identity: {round(percent_identity, 2)}%")
    print(f"This is the alignment score: {filled_matrix[-1][-1]}\n\n")
    return filled_matrix


def main():
    seq1 = "THISLINE"
    seq2 = "ISALIGNED"
    # seq3: GPA1_ARATH
    seq3 = ("MGLLCSRSRHHTEDTDENTQAAEIERRIEQEAKAEKHIRKLLLLGAGESGKSTIFKQIKLLFQ"
            "TGFDEGELKSYVPVIHANVYQTIKLLHDGTKEFAQNETDSAKYMLSSESIAIGEKLSEIGGRLDYPRLTKD"
            "IAEGIETLWKDPAIQETCARGNELQVPDCTKYLMENLKRLSDINYIPTKEDVLYARVRTTGVVEIQFSPVG"
            "ENKKSGEVYRLFDVGGQRNERRKWIHLFEGVTAVIFCAAISEYDQTLFEDEQKNRMMETKELFDWVLKQPC"
            "FEKTSFMLFLNKFDIFEKKVLDVPLNVCEWFRDYQPVSSGKQEIEHAYEFVKKKFEELYYQNTAPDRVDRV"
            "FKIYRTTALDQKLVKKTFKLVDETLRRRNLLEA")
    # seq4: GPA1 BRANA
    seq4 = ( \
        "MGLLCSRSRHHTEDTDENAQAAEIERRIEQEAKAEKHIRKLLLLGAGESGKSTIFKQASS"
        "DKRKIIKLLFQTGFDEGELKSYVPVIHANVYQTIKLLHDGTKEFAQNETDPAKYTLSSEN"
        "MAIGEKLSEIGARLDYPRLTKDLAEGIETLWNDPAIQETCSRGNELQVPDCTKYLMENLK"
        "RLSDVNYIPTKEDVLYARVRTTGVVEIQFSPVGENKKSGEVYRLFDVGGQRNERRKWIHL"
        "FEGVTAVIFCAAISEYDQTLFEDEQKNRMMETKELFDWVLKQPCFEKTSIMLFLNKFDIF"
        "EKKVLDVPLNVCEWFRDYQPVSSGKQEIEHAYEFVKKKFEELYYQNTAPDRVDRVFKIYR"
        "TTALDQKLVKKTFKLVDETLRRRNLLEAGLL")

    matrix = create_matrix(elongate_sequences(seq1), elongate_sequences(seq2))
    fill_matrix(elongate_sequences(seq1), elongate_sequences(seq2), matrix, -8, -8) # Fig 5.9
    fill_matrix(elongate_sequences(seq1), elongate_sequences(seq2), matrix, -4, -4) # Fig 5.11
    fill_matrix(elongate_sequences(seq1), elongate_sequences(seq2), matrix, 0, -8) # Fig 5.12

    matrix_2 = create_matrix(elongate_sequences(seq3), elongate_sequences(seq4))
    fill_matrix(elongate_sequences(seq3), elongate_sequences(seq4), matrix_2, -1, -5)
    fill_matrix(elongate_sequences(seq3), elongate_sequences(seq4), matrix_2, -10, -5)


if __name__ == "__main__":
    main()
