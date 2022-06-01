#!/usr/bin/env python3

"""
Author: Noël Huang

Description: this is a script to ...
"""
# Import statements
import sys
import random

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

def sample(events):
    """Return a key from dict based on the probabilities

    events: dict of {key: probability}, probabilities can also be weights.
    """
    
    pick = random.choices(list(events.keys()),list(events.values()))[0]
    return pick


def parse(file_name):
    """
    Parses a FASTA file containing queries.

    :param file_name: Name of the query FASTA file.
    :return: (list of strings) List of queries
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


def main():
    seq_dict = parse("test_large.fasta")
    for key, value in seq_dict.items():
        print(f"{key}\n{value}")



if __name__ == "__main__":

    # implement main code here
    infile = 'test.fasta'

    main()

    # Put function calls, print statements etc. to answer the questions here
    # When we run your script we should see the answers on screen (or file)
