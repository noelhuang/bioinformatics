#!/usr/bin/env python3

"""
Author: Noël Huang

Description: this is a script to build and train an HMM model based on
previously aligned amino acid sequences.
The trained HMM model can be used to randomly generate an amino acid sequence,
which is similar to the sequences that the model was trained on.

Usage: CMD line: py HMM_model.py <test.fasta> <test_large.fasta>
Where <test.fasta> = fasta file with test amino acid sequences.
<test_large.fasta> = fasta file with large amino acid sequences.
"""
# Import statements
import sys
import random

# Function definitions

# Background amino acid probabilities
pa = {'A': 0.074, 'C': 0.025, 'D': 0.054, 'E': 0.054, 'F': 0.047, 'G': 0.074, \
      'H': 0.026, 'I': 0.068, 'L': 0.099, 'K': 0.058, 'M': 0.025, 'N': 0.045, \
      'P': 0.039, 'Q': 0.034, 'R': 0.052, 'S': 0.057, 'T': 0.051, 'V': 0.073, \
      'W': 0.013, 'Y': 0.034}


class HMM():
    """HMM object to store an HMM model

    This object is designed to keep track of all HMM states, emissions, and
    transitions. It may be used in your implementation, but may also be
    ignored, and replaced by a data structure of choice
    """
    # Emission probabilities for the match and insert states
    e_m = [];
    e_i = pa;

    # Transition probabilities from/to matches, inserts and deletions
    t_mm = [];
    t_mi = [];
    t_md = [];
    t_im = [];
    t_ii = [];
    t_id = [];
    t_dm = [];
    t_di = [];
    t_dd = [];

    def __init__(self, nmatches):
        """Initialize HMM object with number of match states

        nmatches: int, number of match states
        """

        self.nmatches = nmatches

        self.e_m = [dict(pa) for i in range(0, nmatches)]
        for i in range(0, nmatches):
            for j in pa.keys():
                self.e_m[i][j] = 0.0
        self.e_i = pa;

        self.t_mm = [0.0 for i in range(0, nmatches + 1)]
        self.t_mi = [0.0 for i in range(0, nmatches + 1)]
        self.t_md = [0.0 for i in range(0, nmatches + 1)]
        self.t_im = [0.0 for i in range(0, nmatches + 1)]
        self.t_ii = [0.0 for i in range(0, nmatches + 1)]
        self.t_id = [0.0 for i in range(0, nmatches + 1)]
        self.t_dm = [0.0 for i in range(0, nmatches + 1)]
        self.t_di = [0.0 for i in range(0, nmatches + 1)]
        self.t_dd = [0.0 for i in range(0, nmatches + 1)]

    def print_probability(self):
        matrix = [self.t_mm, self.t_mi, self.t_md,
         self.t_im, self.t_ii, self.t_id,
         self.t_dm, self.t_di, self.t_dd]
        # print(matrix)
        return matrix

def sample(events):
    """Return a key from dict based on the probabilities

    events: dict of {key: probability}, probabilities can also be weights.
    """

    pick = random.choices(list(events.keys()), list(events.values()))[0]
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
    Determines which indices in a multiple sequence alignment are match state

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

    # print(match_states)
    # print(amino_count_list)
    return match_states


def reduce_alignment(match_states, seq_dict):
    """
    Using the positions at which a Match is known, produces a reduced alignment

    :param match_states: (list) list whose length is equal to sequence length.
    At each position in the list, presence of a match is indicated by an 'M',
    absence of match is indicated by None.
    :param seq_dict: (dict), with key = str = sequence name,
                            value = str = amino acid sequence.
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
    return reduced_alignments


def count_transitions(reduced_alignments, match_states, seq_dict):
    """
    Counts how many times each transition has occurred for each match position.

    :param reduced_alignments: (dict) reduced alignments, where
                                key = (str), name of sequence
                                value = (str), reduced sequence
    :param match_states: list whose length is equal to sequence length.
    At each position in the list, presence of a match is indicated by an 'M',
    absence of match is indicated by None.
    :param seq_dict: (dict), with key = str = sequence name,
                            value = str = amino acid sequence.
    :return: transition_matrix:
            model, HMM model object with counts of transitions for each
            position in the reduced alignment
    """
    nmatches = len(list(reduced_alignments.values())[0])
    model = HMM(nmatches)

    transition_matrix = [model.t_mm, model.t_mi, model.t_md,
                         model.t_im, model.t_ii, model.t_id,
                         model.t_dm, model.t_di, model.t_dd]

    match_states.append('M')

    sequences = list(seq_dict.values())
    for sequence in sequences:
        sequence = sequence + 'E'  # Add an E at the end to signify end column
        pos_tracker = 0
        state_tracker = 'M'
        for i in range(0, len(sequence)):
            if state_tracker == 'M':  # Start from match col
                if match_states[i] == 'M':  # If going into a match col
                    if sequence[i] == '-':  # If match col contains gap
                        model.t_md[pos_tracker] += 1
                        pos_tracker += 1
                        state_tracker = 'D'
                    elif sequence[i] != '-':  # If match col contains aa
                        model.t_mm[pos_tracker] += 1
                        pos_tracker += 1
                        state_tracker = 'M'
                elif match_states[i] is None:  # If going into insert col
                    # if sequence[i] == '-':
                    #     print("Nothing happens")
                    if sequence[i] != '-':
                        model.t_mi[pos_tracker] += 1
                        state_tracker = None
            elif state_tracker is None:  # Start from insertion column
                if match_states[i] == 'M':  # If going into match col
                    if sequence[i] == '-':  # If deletion found
                        model.t_id[pos_tracker] += 1
                        pos_tracker += 1
                        state_tracker = 'D'
                    elif sequence[i] != '-':  # If aa found
                        model.t_im[pos_tracker] += 1
                        pos_tracker += 1
                        state_tracker = 'M'
                elif match_states[i] is None:  # If going into insertion col
                    # if sequence[i] == '-':
                    #     print("Nothing happens")
                    if sequence[i] != '-':
                        model.t_ii[pos_tracker] += 1
            elif state_tracker == 'D':  # If starting from deletion state
                if match_states[i] == 'M':  # And if going into match column
                    if sequence[i] == '-':  # If gap is found (aka deletion)
                        model.t_dd[pos_tracker] += 1
                        pos_tracker += 1
                    elif sequence[i] != '-':  # If amino acid is found
                        model.t_dm[pos_tracker] += 1
                        pos_tracker += 1
                        state_tracker = 'M'
                elif match_states[i] is None:  # And if going into insertion col
                    # if sequence[i] == '-':
                    #     print("Nothing happens")
                    if sequence[i] != '-':
                        model.t_di[pos_tracker] += 1
                        state_tracker = None

    # print(sequences)
    # for transition in transition_matrix:
    #     print(transition)
    # print(model.e_m)
    return transition_matrix, model


def calculate_probabilities(transition_matrix, model):
    """
    Calculates probabilities for every transition in the HMM model

    :param transition_matrix: list of list of floats. The list of floats
    represent a row in the HMM model, there are 9 rows, which represent (from
    1 to 9) the following transitions; mm, mi, md, im, ii, id, dm, di, dd.
    The column number in the rows represents the match state index.
    :param model: (object) HMM model, created with the HMM class.
    :return: (object) HMM model with calculated transmission probabilities
    """
    probability_matrix = [x[:] for x in transition_matrix]
    out_degree_matrix = [[], [], []]
    for j in range(0, len(transition_matrix[0])):
        out_degree_m = 0
        out_degree_i = 0
        out_degree_d = 0
        for i in range(0, 3):
            out_degree_m += transition_matrix[i][j]
        out_degree_matrix[0].append(out_degree_m)

        for i in range(3, 6):
            out_degree_i += transition_matrix[i][j]
        out_degree_matrix[1].append(out_degree_i)

        for i in range(6, 9):
            out_degree_d += transition_matrix[i][j]
        out_degree_matrix[2].append(out_degree_d)
    # print(transition_matrix)
    for j in range(0, len(probability_matrix[0])):
        for i in range(0, 3):
            if out_degree_matrix[0][j] != 0:
                probability_matrix[i][j] = probability_matrix[i][j] \
                                           / out_degree_matrix[0][j]
        for i in range(3, 6):
            if out_degree_matrix[1][j] != 0:
                probability_matrix[i][j] = probability_matrix[i][j] \
                                           / out_degree_matrix[1][j]
        for i in range(6, 9):
            if out_degree_matrix[2][j] != 0:
                probability_matrix[i][j] = probability_matrix[i][j] \
                                           / out_degree_matrix[2][j]

    for j in range(0, len(model.t_mm)):
        model.t_mm[j] = probability_matrix[0][j]
        model.t_mi[j] = probability_matrix[1][j]
        model.t_md[j] = probability_matrix[2][j]
        model.t_im[j] = probability_matrix[3][j]
        model.t_ii[j] = probability_matrix[4][j]
        model.t_id[j] = probability_matrix[5][j]
        model.t_dm[j] = probability_matrix[6][j]
        model.t_di[j] = probability_matrix[7][j]
        model.t_dd[j] = probability_matrix[8][j]
    return model


def calculate_emission_probability(reduced_alignments, model):
    """
    Calculates emission probabilities in an HMM model

    :param reduced_alignments: (dict) reduced alignments, where
                                key = (str), name of sequence
                                value = (str), reduced sequence
    :param model: (object) HMM model containing not-yet-calculated emission
                            probabilities
    :return: (object) HMM model containing calculated emission probabilities

    HMM model was created using HMM class
    """
    amino_counter = 0
    amino_count_list = []
    reduced_alignments = list(reduced_alignments.values())
    alignment_length = len(reduced_alignments[0])
    for i in range(0, alignment_length):
        for sequence in reduced_alignments:
            if sequence[i] != "-":
                amino_counter = amino_counter + 1
        amino_count_list.append(amino_counter)
        amino_counter = 0

    for sequence in reduced_alignments:
        for i, letter in enumerate(sequence):
            if letter != '-':
                model.e_m[i][letter] += 1

    for pos, dictionary in enumerate(model.e_m):
        for key, value in dictionary.items():
            dictionary[key] = value / amino_count_list[pos]
    # print(model.e_m)
    return model


def generate_sequence(model, reduced_alignments):
    """
    Randomly generates an amino acid sequence using an HMM model.

    :param model: (object) HMM model containing transmission and emission
                            probabilities.
    :param reduced_alignments: (dict) reduced alignments, where
                                key = (str), name of sequence
                                value = (str), reduced sequence

    :return: (str) a sequence that was randomly generated using the HMM model.
    """
    nmatches = len(list(reduced_alignments.values())[0])
    transition_dict = {'mm': 0.0, 'mi': 0.0, 'md': 0.0,
                       'im': 0.0, 'ii': 0.0, 'id': 0.0,
                       'dm': 0.0, 'di': 0.0, 'dd': 0.0}
    transition_prob_dicts = [dict(transition_dict) for i in range(0, nmatches+1)]
    for i in range(0, len(transition_prob_dicts)):
        transition_prob_dicts[i]['mm'] = model.t_mm[i]
        transition_prob_dicts[i]['mi'] = model.t_mi[i]
        transition_prob_dicts[i]['md'] = model.t_md[i]
        transition_prob_dicts[i]['im'] = model.t_im[i]
        transition_prob_dicts[i]['ii'] = model.t_ii[i]
        transition_prob_dicts[i]['id'] = model.t_id[i]
        transition_prob_dicts[i]['dm'] = model.t_dm[i]
        transition_prob_dicts[i]['di'] = model.t_di[i]
        transition_prob_dicts[i]['dd'] = model.t_dd[i]

    # print(len(transition_prob_dicts), transition_prob_dicts)
    #
    # print("model e_m",len(model.e_m), model.e_m)
    sequence = ''
    pos_tracker = 0
    state_tracker = 'M'
    match_keys = ['mm', 'mi', 'md']
    insertion_keys = ['im', 'ii', 'id']
    deletion_keys = ['dm', 'di', 'dd']

    while pos_tracker <= nmatches:
        # print('entering with pos tracker:', pos_tracker)
        if state_tracker == 'M':
            pick = random.choices(match_keys,
                                  [transition_prob_dicts[pos_tracker][key]
                                   for key in match_keys])[0]
            if pick == 'mm':
                amino = sample(model.e_m[pos_tracker-1])
                sequence += amino
                pos_tracker += 1
            elif pick == 'mi':
                amino = sample(model.e_i)
                sequence += amino
                state_tracker = None
            elif pick == 'md':
                sequence += '-'
                pos_tracker += 1
                state_tracker = 'D'

        elif state_tracker is None:
            pick = random.choices(insertion_keys,
                                  [transition_prob_dicts[pos_tracker][key]
                                   for key in insertion_keys])[0]
            if pick == 'im':
                amino = sample(model.e_m[pos_tracker-1])
                sequence += amino
                pos_tracker += 1
                state_tracker = 'M'
            elif pick == 'ii':
                amino = sample(model.e_i)
                sequence += amino
            elif pick == 'id':
                sequence += '-'
                pos_tracker += 1
                state_tracker = 'D'

        elif state_tracker == 'D':
            pick = random.choices(deletion_keys,
                                  [transition_prob_dicts[pos_tracker][key]
                                   for key in deletion_keys])[0]
            if pick == 'dm':
                amino = sample(model.e_m[pos_tracker-1])
                sequence += amino
                pos_tracker += 1
                state_tracker = 'M'
            elif pick == 'di':
                amino = sample(model.e_i)
                sequence += amino
                state_tracker = None
            elif pick == 'dd':
                sequence += '-'
                pos_tracker += 1

    # print('helptest', sequence)
    # print('postracker', pos_tracker)
    return sequence


def main():
    ### Question 1
    seq_dict = parse(sys.argv[1])
    match_states = calc_match_states(seq_dict)
    reduced_alignments = reduce_alignment(match_states, seq_dict)
    print("Question 1:")
    print("Number of match states:", len(list(reduced_alignments.values())[0]))

    ### Question 2
    f = open('reduced_fasta.txt', 'w')
    for key, value in reduced_alignments.items():
        f.write(key)
        f.write('\n')
        f.write(value)
        f.write('\n')
    f.close()

    ### Question 3
    transition_matrix, model = \
        count_transitions(reduced_alignments, match_states, seq_dict)
    calculate_probabilities(transition_matrix, model)

    calculate_emission_probability(reduced_alignments, model)
    print('Question 3:')
    print(f" t_mm {model.t_mm}\n t_mi {model.t_mi}\n t_md {model.t_md}\n "
          f"t_im {model.t_im}\n t_ii {model.t_ii}\n t_id {model.t_id}\n "
          f"t_dm {model.t_dm}\n t_di {model.t_di}\n t_dd {model.t_dd}\n ")
    print(f'Match emission probability:\n{model.e_m}')

    ### Question 4
    print("Question 4:")
    for i in range(0, 10):
        sequence = generate_sequence(model, reduced_alignments)
        print(f"Seq. no. {i+1}:", sequence)

    ### Question 5-1
    seq_dict = parse(sys.argv[2])
    match_states = calc_match_states(seq_dict)
    reduced_alignments = reduce_alignment(match_states, seq_dict)

    ### Question 5-2
    f = open('reduced_fasta_question_5_2.txt', 'w')
    for key, value in reduced_alignments.items():
        f.write(key)
        f.write('\n')
        f.write(value)
        f.write('\n')
    f.close()

    ### Question 5-3
    transition_matrix, model = \
        count_transitions(reduced_alignments, match_states, seq_dict)
    calculate_probabilities(transition_matrix, model)

    calculate_emission_probability(reduced_alignments, model)
    print('Question 5-3:')

    print("First transmissions")
    print(f" t_mm {model.t_mm[0]}\n t_mi {model.t_mi[0]}\n t_md {model.t_md[0]}\n "
          f"t_im {model.t_im[0]}\n t_ii {model.t_ii[0]}\n t_id {model.t_id[0]}\n "
          f"t_dm {model.t_dm[0]}\n t_di {model.t_di[0]}\n t_dd {model.t_dd[0]}\n ")

    print("Last transmissions")
    print(f" t_mm {model.t_mm[-1]}\n t_mi {model.t_mi[-1]}\n t_md {model.t_md[-1]}\n "
          f"t_im {model.t_im[-1]}\n t_ii {model.t_ii[-1]}\n t_id {model.t_id[-1]}\n "
          f"t_dm {model.t_dm[-1]}\n t_di {model.t_di[-1]}\n t_dd {model.t_dd[-1]}\n ")

    print(f'First match emission probability:\n{model.e_m[0]}')
    print(f'Last match emission probability:\n{model.e_m[-1]}')


    ### Question 5-4
    print("Question 4:")
    for i in range(0, 10):
        sequence = generate_sequence(model, reduced_alignments)
        print(f"Seq. no. {i+1}:", sequence)


if __name__ == "__main__":
    main()

