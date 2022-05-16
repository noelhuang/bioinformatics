"""
Name: Noël Huang
Reg. nr: 1032641
Python 3 script to construct a de Bruijn graph from DNA string inputs.

Usage: CMD line: py BIF007.py <edge list format file> <number of starting vertex>
Args:
    <dna strings>: A collection of up to 1000 (possibly repeating) DNA strings of equal length
    (length =< 50 bp) corresponding to a set S of (k+1)-mers.
"""

import sys


def parse(file_name):
    """
    Parses a file and returns data in list format, removes duplicates and also makes reverse complement set.

    Args:
        file_name: Name of the file which contains a collection of up to 1000 (possibly repeating)
        DNA strings of equal length
        (length =< 50 bp) corresponding to a set S of (k+1)-mers.

    Returns: List format of all DNA strings of equal length in the file without duplicates.

    """
    with open(file_name) as f:
        parsed_data = f.read()
        parsed_data = parsed_data.split("\n")
        if parsed_data[-1] == "":
            del parsed_data[-1]

        parsed_data = set(parsed_data)
        parsed_data = list(parsed_data)
    return parsed_data


def make_rc(a_string):
    rc_dict = {"A": "T", "C": "G", "G": "C", "T": "A"}
    complementary_list = []
    for i in range(0, len(a_string)):
        complementary_list.append(rc_dict[a_string[i]])
    reverse_complementary = [i for i in complementary_list[-1::-1]]
    joined = "".join(reverse_complementary)
    print(joined)
    return joined


def add_rc_to_list(parsed_data):
    for i in range(0, len(parsed_data)):
        parsed_data.append(make_rc(parsed_data[i]))
    parsed_data = set(parsed_data)
    parsed_data = list(parsed_data)
    return parsed_data


def make_edges(complete_set):
    edge_list = []
    for dna_sequence in complete_set:
        node_u = dna_sequence[0:-1]
        node_v = dna_sequence[1:]
        edge = [node_u, node_v]
        edge_list.append(edge)
    edge_list.sort()
    print(edge_list)
    return edge_list


def main():
    parsed_data = parse(sys.argv[1])
    complete_set = add_rc_to_list(parsed_data)
    edges = make_edges(complete_set)

    neat_values = "\n".join(str(x) for x in edges)
    print(neat_values)
    neat_values = neat_values.replace("'", "")
    neat_values = neat_values.replace("[", "(")
    neat_values = neat_values.replace("]", ")")
    print(type(neat_values))
    print(neat_values)

    with open("output_bruijn.txt", "w") as f:
        f.write(neat_values)


if __name__ == "__main__":
    main()

