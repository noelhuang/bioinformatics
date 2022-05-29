#!/usr/bin/env python3
"""
Author: Noël Huang
Student number: 1032641

Description:Implementation of the SSAHA algorithm.

Usage: CMD line:  py SSAHA.py <arabidopsis genome FASTA> <query FASTA>
Args: <arabidopsis genome FASTA>: FASTA file of the arabidopsis genome.
<query FASTA>: FASTA file of the queries you want to search for.
"""

import sys


def construct_hash_table(sequences, k):
    """
    Constructs a hash table given a database of sequences and word length k.

    :param sequences: list of sequences which form the database which is going
    to be searched in.
    :param k: (int) length of k-tuples.
    :return: (dict) hash_table. A dictionary which contains all
    non-overlapping k-tuples that occur in the database as keys. The values
    are a list of tuples. Each tuple (i, j) represents the index of
    the sequence within the database, while the j represents the index of the
    beginning of the k-tuple within sequence i.

    Example: Hash table: {'AA': [(1, 18)]} indicates that the k-tuple AA
    occurs in sequence 1, at base 18 in that sequence.
    """
    hash_table = {}
    for i in range(0, len(sequences)):
        for j in range(0, len(sequences[i]) - k + 1, k):
            ktuple_in_sequence = sequences[i][j:j + k]
            if ktuple_in_sequence in hash_table:
                hash_table[ktuple_in_sequence].append((i, j))
            else:
                hash_table[ktuple_in_sequence] = [(i, j)]
    hash_table = {key: value for key, value in sorted(hash_table.items())}
    return hash_table


def sequence_search(query, hash_table, k):
    """
    Finds all hits between query and database, returns sorted list of all hits.

    :param query: (str) the DNA sequence you want to match to the database.
    :param hash_table: A dictionary which contains all
    non-overlapping k-tuples that occur in the database as keys. The values
    are a list of tuples. Each tuple (i, j) represents the index of
    the sequence within the database, while the j represents the index of the
    beginning of the k-tuple within sequence i.
    :param k: (int) length of k-tuples.
    :return: A sorted list of tuples. Each tuple represents a match (hit)
    between query and database. Tuple consists of three integers (i, j, k)
    where:
    i = index (index of sequence in database)
    j = shift (number of bases that the query had to be shifted to align with
    the database sequence)
    k = offset (index of first base in the k-tuple in the database sequence)
    """
    hits = []
    for t in range(0, len(query) - k + 1):
        query_ktuple = query[t:t + k]
        if query_ktuple in hash_table:
            positions = hash_table[query_ktuple]
            for position in positions:
                hit = (position[0], position[1] - t, position[1])
                hits.append(hit)
    hits.sort()
    master_list = hits
    print("First and last hits:", master_list[0], master_list[-1])
    return master_list


def find_longest_match(master_list):
    """
    Finds the longest run of hits in the sorted hit list.

    :param master_list: A sorted list of tuples. Each tuple represents a match (hit)
    between query and database. Tuple consists of three integers (i, j, k)
    where:
    i = index (index of sequence in database)
    j = shift (number of bases that the query had to be shifted to align with
    the database sequence)
    k = offset (index of first base in the k-tuple in the database sequence)
    :return: (list of tuples) cons_hits, a list of tuples which
    represent the longest run of consecutive hits
    """
    if len(master_list) > 0:
        max_consecutive_count = 1
    else:
        max_consecutive_count = 0
    consecutive_count = 1
    indices_of_alignments = []
    start_index_align = None

    for n in range(0, len(master_list) - 1):
        if master_list[n][0:2] == master_list[n + 1][0:2]:
            consecutive_count = consecutive_count + 1
            if consecutive_count > max_consecutive_count:
                max_consecutive_count = consecutive_count
                start_index_align = n - consecutive_count + 2
        else:
            indices_of_alignments.append((n - consecutive_count + 1, n))
            consecutive_count = 1

    if max_consecutive_count > 1:
        print("Total number of hits:", len(master_list))
        print("Maximum consecutive hits:", max_consecutive_count)
        end_index_align = start_index_align + max_consecutive_count - 1
        print("Start, end index:", (start_index_align, end_index_align))
        cons_hits = master_list[start_index_align:end_index_align + 1]
        return cons_hits
    else:
        print("Longest alignment was only 1 k-tuple long")
        cons_hits = master_list[0]
        print("Maximum number of consecutive hits was 1:\n"
              "Index, shift, offset of the first hit", cons_hits)
        return [cons_hits]


def print_alignment(query, sequences, cons_hits):
    """
    Gives the longest exact alignment between a query and a database.

    :param query: (str) DNA string that you want to match to database.
    :param sequences: (list of strings) database, contains DNA sequence strings
    :param cons_hits: (list of tuples) a list of tuples which
    represent the longest run of consecutive hits
    :return: (list of strings) list containing 3 strings.
    First string: DNA sequence in database which was matched with
    Second string: alignment indicator string
    Third string: Query DNA sequence which was searched for in database.
    """
    print("Index, shift, offset of consecutive hit(s):", cons_hits)
    index = cons_hits[0][0]
    initial_offset = cons_hits[0][2]
    shift = cons_hits[0][1]
    correct_sequence = sequences[index]
    prefix_length = len(query) // 2
    suffix_length = len(query) // 2
    abbreviation = correct_sequence[initial_offset
                                    - prefix_length:initial_offset
                                                    + len(query)
                                                    + suffix_length]
    database_dna = f"{abbreviation}"
    base_number_t = initial_offset - shift
    query_aligned = ((prefix_length - base_number_t) * " ") + query
    while len(query_aligned) < len(database_dna):
        query_aligned = query_aligned + " "

    alignment_indicator = ""
    for i in range(0, len(query_aligned)):
        if query_aligned[i] == " ":
            alignment_indicator += " "
        elif query_aligned[i] == database_dna[i]:
            alignment_indicator += "|"
        else:
            alignment_indicator += "."

    alignment = [database_dna, alignment_indicator, query_aligned]
    print("This is the alignment:")

    line_length = 79  # This block just makes prettier print output
    lines_0 = [database_dna[i:i + line_length] +
               '\n' for i in range(0, len(database_dna), line_length)]
    lines_1 = [alignment_indicator[i:i + line_length] +
               '\n' for i in range(0, len(alignment_indicator), line_length)]
    lines_2 = [query_aligned[i:i + line_length] +
               '\n' for i in range(0, len(query_aligned), line_length)]

    for i in range(0, len(lines_0)):
        print(f"{lines_0[i]}{lines_1[i]}{lines_2[i]}")

    print("\n")
    return alignment


def parse_genome(file_name):
    """
    Parses a FASTA genome file.

    :param file_name: (str) name of the genome FASTA file.
    :return: (str) genome DNA sequence
    """
    with open(file_name) as file:
        data = file.read()
        data = data.split("\n")
        sequence_count = 0
        for line in data:
            if ">" in line:
                sequence_count += 1
                data.remove(line)
            if line == '':
                data.remove(line)
            elif len(line) == 0:
                data.remove(line)

        data = ''.join(data)
    print("Combined length of sequences (base pairs):", len(data))
    print("Sequence count", sequence_count)
    return [data]


def parse_query(file_name):
    """
    Parses a FASTA file containing queries.

    :param file_name: Name of the query FASTA file.
    :return: (list of strings) List of queries
    """
    with open(file_name) as file:
        data = file.read()
        data = data.split("\n")
        for line in data:
            if ">" in line:
                query_name = line
                data.remove(line)

        queries = [""]

        for i in range(0, len(data)):
            if len(data[i]) == 0:
                queries.append("")
            else:
                queries[-1] = queries[-1] + data[i]
        return queries


def main():
    sequences = ["GTGACGTCACTCTGAGGATCCCCTGGGTGTGG",
                 "GTCAACTGCAACATGAGGAACATCGACAGGCCCAAGGTCTTCCT",
                 "GGATCCCCTGTCCTCTCTGTCACATA"]
    query = "TGCAACAT"
    k = 2

    print("--Start Question 1--")
    hash_table = construct_hash_table(sequences, k)
    print("Hash table:", hash_table)
    print("--End Question 1--\n")

    print("--Start Question 2--")
    master_list = sequence_search(query, hash_table, k)
    cons_hits = find_longest_match(master_list)
    print("Sorted hit list:", master_list)
    print("--End Question 2--\n")

    print("--Start Question 3--")
    print_alignment(query, sequences, cons_hits)
    print("--End Question 3--\n")

    print("--Start Question 4--")
    genome_file = sys.argv[1]  # File name for arabidopsis genome
    arabidopsis_data = parse_genome(genome_file)
    print("--End Question 4--\n")

    print("--Start Question 5--")
    k = 13
    hash_table_arabidopsis = construct_hash_table(arabidopsis_data, k)
    print("Length hash table:", len(hash_table_arabidopsis))
    print("--End Question 5--\n")

    print("--Start Question 6--")
    query_file = sys.argv[2]
    queries = parse_query(query_file)
    for query in queries:
        print("Displaying results for this query:", query)
        master_list = sequence_search(query, hash_table_arabidopsis, k)
        cons_hits = find_longest_match(master_list)
        print_alignment(query, arabidopsis_data, cons_hits)
    print("--End Question 6--\n")


if __name__ == "__main__":
    main()
