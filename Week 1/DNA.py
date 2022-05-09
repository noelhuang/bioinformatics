# Student name: Noël Huang
# Student number: 1032641

import sys

file_name = sys.argv[1]


def parse_file():
    with open(file_name) as data:
        sequence = data.read()
        print("This is the input sequence: \n" + f"{sequence}")
        return sequence


def gc_content(sequence):
    gc_count = 0
    gc_list = ["g", "c", "G", "C"]
    for i in sequence:
        if i in gc_list:
            gc_count = gc_count + 1
    gc_content = gc_count / len(sequence)
    return gc_content


def reverse_complement(sequence):
    base_num = {'g': '1',
                't': '2',
                'a': '3',
                'c': '4',
                'G': '5',
                'T': '6',
                'A': '7',
                'C': '8'}

    num_base = {'1': 'c',
                '2': 'a',
                '3': 't',
                '4': 'g',
                '5': 'C',
                '6': 'A',
                '7': 'T',
                '8': 'G'}
    for key, value in base_num.items():
        sequence = sequence.replace(key, value)
    for key, value in num_base.items():
        sequence = sequence.replace(key, value)

    return sequence


def main():
    sequence = parse_file()
    print(f"Input length is: {len(sequence)} bp")
    print(f"GC content is: {round(gc_content(sequence) * 100, 2)}%")
    print(f"Reverse complement is: \n{reverse_complement(sequence)}")


if __name__ == "__main__":
    main()
