"""
Name: Noël Huang
Reg. nr: 1032641
Python 3 script to count number of A C G and T bases in a sequence (string).

Usage: CMD line: py BIF002.py <DNA sequence file>

Args: <DNA sequence file>: Name of a .txt file containing uppercase DNA sequence.
"""
import sys


def parse(file_name):
    """
    Parses a file and returns contents as a string

    :param file_name: (str) name of the file to parse
    :return: data: (str) contents of the file
    """
    with open(file_name) as data:
        data = data.read()
    return data


def count_bases(some_string):
    """
    Counts the number of ACGT bases in a string

    :param some_string: (str) String containing a sequence of bases ACGT (capital letters)
    :return: A string of four numbers separated by a string, which represent number of A, C, G,
    and T bases respectively.
    """
    # List of base count (order: A,C,G,T)
    base_dic = {"A": 0, "C": 0, "G": 0, "T": 0}
    for i in range(0, len(some_string)):
        if some_string[i] in base_dic:
            base_dic[some_string[i]] = base_dic[some_string[i]] + 1
    return base_dic


def main():
    """
    Prints A C G T count separated by space

    """
    file_name = sys.argv[1]
    sequence = parse(file_name)
    base_dic = count_bases(sequence)
    print(base_dic['A'], base_dic['C'], base_dic['G'], base_dic['T'])


if __name__ == "__main__":
    main()
