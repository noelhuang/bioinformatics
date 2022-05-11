"""
Name: Noël Huang
Reg. nr: 1032641
Python 3 script to count number of 4-mers in a sequence (string).
Prints an array of integers, which represent the number of times the alphabetically sorted 4-mers appear in the string.

Usage: CMD line: py BIF003.py <DNA sequence file>
Args: <DNA sequence file>: Name of a .FAS or .FASTA file containing uppercase DNA sequence. .txt is also allowed
"""
import sys


def generate_4mer():
    """
    Produces a list of all possible 4-mers using bases A, C, G, and T, sorted lexicographically.

    :return kmer4321: A list of all possible 4-mers using bases A, C, G, and T, sorted lexicographically.
    """
    base_dic = {1: "A", 2: "C", 3: "G", 4: "T"}
    kmer1 = []
    for i in range(1, 5):
        base1 = base_dic[i]
        kmer1.append(base1)
    kmer21 = []
    for i in range(1, 5):
        base2 = base_dic[i]
        for base1 in kmer1:
            base21 = base2 + base1
            kmer21.append(base21)
    kmer321 = []
    for i in range(1, 5):
        base3 = base_dic[i]
        for base21 in kmer21:
            base321 = base3 + base21
            kmer321.append(base321)
    kmer4321 = []
    for i in range(1, 5):
        base4 = base_dic[i]
        for base321 in kmer321:
            base4321 = base4 + base321
            kmer4321.append(base4321)
    print(kmer4321)
    return kmer4321


def parse(file_name):
    """
    Parses a FASTA file, removes annotation (">annotation") and returns remaining contents as a string

    :param file_name: (str) name of the file to parse, must be a .FAS or .FASTA file
    :return: data: (str) contents of the file as string. In this case string of DNA sequence
    """
    with open(file_name) as data:
        data = data.read()
        data = data.split("\n")
        for element in data:
            if ">" in element:
                data.remove(element)
    data = "".join([item for item in data])
    return data


def find_4mers_in_dna(sequence):
    """
    Finds and returns all 4-mers in a given DNA sequence, in order from left to right.

    :param sequence: String of DNA sequence with bases A C G T
    :return: A list of all the 4-mers in the DNA sequence, in order from left to right.
    """
    fourmers_in_dna = []
    for i in range(0, len(sequence)-3):
        fourmer = sequence[i] + sequence[i+1] + sequence[i+2] + sequence[i+3]
        fourmers_in_dna.append(fourmer)

    print(fourmers_in_dna)
    return fourmers_in_dna


def main():
    file_name = sys.argv[1]
    sequence = parse(file_name)
    kmer4321 = generate_4mer()
    fourmers_in_dna = find_4mers_in_dna(sequence)

    count_list = [0 for item in kmer4321]

    for fourmer in fourmers_in_dna:
        if fourmer in kmer4321:
            index = kmer4321.index(fourmer)
            count_list[index] = count_list[index] + 1

    print(" ".join(str(x) for x in count_list))


if __name__ == "__main__":
    main()
