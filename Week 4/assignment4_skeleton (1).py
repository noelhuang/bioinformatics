#!/usr/bin/env python3

"""
Author: Noël Huang

Description: this is a script to ...
"""
# Import statements
import sys
import random
import numpy as np

# Function definitions

# Background amino acid probabilities
pa = { 'A':0.074, 'C':0.025, 'D':0.054, 'E':0.054, 'F':0.047, 'G':0.074,\
    'H':0.026, 'I':0.068, 'L':0.099, 'K':0.058, 'M':0.025, 'N':0.045,\
    'P':0.039, 'Q':0.034, 'R':0.052, 'S':0.057, 'T':0.051, 'V':0.073,\
    'W':0.013, 'Y':0.034 }


class HMM():
    """HMM object to store an HMM model

    This object is designed to keep track of all HMM states, emissions, and
    transitions. It may be used in your implementation, but may also be
    ignored, and replaced by a data structure of choice
    """
    # Emission probabilities for the match and insert states
    e_m   = []; e_i   = pa;

    # Transition probabilities from/to matches, inserts and deletions
    t_mm = []; t_mi = []; t_md = [];
    t_im = []; t_ii = []; t_id = [];
    t_dm = []; t_di = []; t_dd = [];

    def __init__(self,nmatches):
        """Initialize HMM object with number of match states

        nmatches: int, number of match states
        """

        self.nmatches = nmatches

        self.e_m   = [dict(pa) for i in range(0,nmatches)]
        for i in range(0,nmatches):
            for j in pa.keys():
                self.e_m[i][j] = 0.0
        self.e_i   = pa;

        self.t_mm  = [0.0 for i in range(0,nmatches+1)]
        self.t_mi  = [0.0 for i in range(0,nmatches+1)]
        self.t_md  = [0.0 for i in range(0,nmatches+1)]
        self.t_im  = [0.0 for i in range(0,nmatches+1)]
        self.t_ii  = [0.0 for i in range(0,nmatches+1)]
        self.t_id  = [0.0 for i in range(0,nmatches+1)]
        self.t_dm  = [0.0 for i in range(0,nmatches+1)]
        self.t_di  = [0.0 for i in range(0,nmatches+1)]
        self.t_dd  = [0.0 for i in range(0,nmatches+1)]


# class HMM_pos:
#
#     def __init__(self, sequences):
#         for i in range(0, len(sequences[0])):
#             if i == 0:
#                 self.t_mm =
#                 self.t_mi =
#                 self.t_md =
#                 self.t_im =
#                 self.t_ii =
#                 self.t_id =
#                 self.t_dm =
#                 self.t_di =
#                 self.t_dd =

def sample(events):
    """Return a key from dict based on the probabilities

    events: dict of {key: probability}, probabilities can also be weights.
    """
    
    pick = random.choices(list(events.keys()),list(events.values()))[0]
    return pick


def parse(file_name):
    """
    Parses a FASTA file containing DNA/amino acid sequences.

    :param file_name: Name of the FASTA file to be parsed.
    :return: (dict), with key = str = sequence name, value = str = DNA/amino acid sequence.
    """
    seq_dict = {}
    with open(file_name) as file:
        data = file.read()
        data = data.split("\n")
        for line in data:
            if ">" in line:
                seq_name = line
                seq_dict[seq_name] = ""
            if ">" not in line and len(line) > 0:
                seq_dict[seq_name] = seq_dict[seq_name] + line

    return seq_dict


def calc_match_states(seq_dict):
    """
    Determines which positions in a multiple sequence alignment are match states

    :param seq_dict: (dict), with key = str = sequence name,
                            value = str = DNA/amino acid sequence.
    :return: list whose length is equal to sequence length. At each position
    in the list, presence of a match is indicated by an 'M',
    absence of match is indicated by None.

    A position is set to match, if at this position at least half of the
    sequences have an amino acid.
    """
    amino_counter = 0
    amino_count_list = []
    sequences = list(seq_dict.values())

    sequence_len = len(sequences[0])

    number_of_sequences = len(sequences)
    for i in range(0, sequence_len):
        for sequence in sequences:
            if sequence[i] != "-":
                amino_counter = amino_counter + 1
        amino_count_list.append(amino_counter)
        amino_counter = 0

    match_states = [None for i in range(0, sequence_len)]

    for i, count in enumerate(amino_count_list):
        if count >= 0.5 * number_of_sequences:
            match_states[i] = "M"

    print(match_states)
    print(amino_count_list)
    return match_states


def reduce_alignment(match_states, seq_dict):
    """
    Using the positions at which a Match is known, produces a reduced alignment

    :param match_states: (list) list whose length is equal to sequence length.
    At each position in the list, presence of a match is indicated by an 'M',
    absence of match is indicated by None.
    :param seq_dict: (dict), with key = str = sequence name,
                            value = str = DNA/amino acid sequence.
    :return: dictionary with sequence name (str) as keys, and reduced alignment
    sequences as values (str).
    """
    reduced_alignments = {}
    for seq_name, sequence in seq_dict.items():
        for j, char in enumerate(sequence):
            if match_states[j] == 'M':
                if seq_name in reduced_alignments:
                    reduced_alignments[seq_name] += char
                else:
                    reduced_alignments[seq_name] = char
    print(reduced_alignments)
    return reduced_alignments


# def make_match_dict(match_states):
#     match_dict = {}
#     match_indices = []
#     for index, state in enumerate(match_states):
#         if state == 'M':
#             match_indices.append(index)
#     print("match indices", match_indices)
#     for i in range(0, len(match_indices)):
#         if i < len(match_indices)-1:
#             match_dict[match_indices[i]] = match_indices[i+1]
#         if i == len(match_indices)-1:
#             match_dict[match_indices[i]] = None
#
#     return


def count_transitions_emissions(reduced_alignments, match_states, seq_dict):
    nmatches = len(list(reduced_alignments.values())[0])
    model = HMM(nmatches)

    sequences = list(seq_dict.values())
    for sequence in sequences:
        pos_tracker = 0
        state_tracker = match_states[0]
        for pos in range(0, len(sequence)):
            if state_tracker is None: # If last state was insertion
                if sequence[pos] != '-': # And if an amino acid is found at pos
                    state_tracker = match_states[pos] # Update state tracker to insert
                    model.t_



            if current_state == 'I':

    print(sequences)
    print(model.t_mm)




def main():
    seq_dict = parse("develop.fasta")
    for key, value in seq_dict.items():
        print(f"{key}\n{value}")

    match_states = calc_match_states(seq_dict)
    reduced_alignments = reduce_alignment(match_states, seq_dict)

    count_transitions_emissions(reduced_alignments, match_states, seq_dict)


if __name__ == "__main__":

    # implement main code here
    infile = 'develop.fasta'

    main()

    # Put function calls, print statements etc. to answer the questions here
    # When we run your script we should see the answers on screen (or file)
