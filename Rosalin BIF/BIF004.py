"""
Name: Noël Huang
Reg. nr: 1032641
Python 3 script to find all k-mers given an integer k and a DNA sequence (str).

Usage: CMD line: py BIF004.py <DNA sequence file>
Args: <DNA sequence file>: Name of a .txt file containing integer k uppercase and a DNA sequence separated by \n
"""
import sys


def parse(file_name):
    """
    Parses a txt file, returns integer k and DNA sequence.

    :param file_name: (str) name of the file to parse, must be a .txt file with integer k followed by \n DNA sequence
    :return: k, data: integer k denotes k-mer length, data denotes string of DNA sequence.
    """
    with open(file_name) as data:
        data = data.read()
        data = data.split("\n")
        k = int(data[0])
        del data[0]
    data = "".join([item for item in data])
    return k, data


def detect_kmer(k, data):
    """

    :param k: integer k denotes kmer length
    :param data: string data denotes sequence of DNA, capital letters, bases ACGT
    :return: list of k-mers (str), ordered from detection left to right
    """
    kmer_list = []
    for i in range(0, len(data)-k+1):
        kmer_list.append(data[i:i+k])
    return kmer_list


def main():
    file_name = sys.argv[1]
    k, data = parse(file_name)
    kmer_list = detect_kmer(k, data)
    neat_output = "\n".join(str(x) for x in kmer_list)
    print(neat_output)
    f = open("output.txt", "w")
    f.write(neat_output)


if __name__ == "__main__":
    main()
