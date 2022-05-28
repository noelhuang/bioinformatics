import itertools as itertools


# def make_all_ktuples(alphabet, k):
#     ktuples = list(itertools.product(alphabet, repeat=k))
#     for i in range(0, len(ktuples)):
#         ktuples[i] = ''.join(ktuples[i])
#
#     return ktuples


def construct_hash_table(sequences, k):

    # empty_buckets = [None for i in range(len(all_ktuples))]

    hash_table = {}

    for i in range(0, len(sequences)):
        for j in range(0, len(sequences[i])-k+1, k):
            ktuple_in_sequence = sequences[i][j:j+k]
            print("Ktuple in  database sequence:", ktuple_in_sequence)
            if ktuple_in_sequence in hash_table:
                hash_table[ktuple_in_sequence].append((i, j))
            else:
                hash_table[ktuple_in_sequence] = [(i, j)]
    hash_table = {key: value for key, value in sorted(hash_table.items())}
    print(hash_table)
    return hash_table


def sequence_search(query, hash_table, k):
    hits = []
    for t in range(0, len(query)-k+1):
        query_ktuple = query[t:t+k]
        if query_ktuple in hash_table:
            positions = hash_table[query_ktuple]
            for position in positions:
                hit = (position[0], position[1]-t, position[1])
                hits.append(hit)
    hits.sort()
    master_list = hits
    print(master_list)
    return master_list


def find_longest_match(master_list):
    max_consecutive_count = 0
    consecutive_count = 1
    indices_of_alignments = []
    index_longest_align = None

    for n in range(0, len(master_list)-1):
        if master_list[n][0:2] == master_list[n+1][0:2]:
            consecutive_count = consecutive_count + 1
            if consecutive_count > max_consecutive_count:
                max_consecutive_count = consecutive_count
                index_longest_align = n - consecutive_count + 2
        else:
            indices_of_alignments.append((n - consecutive_count + 1, n))
            consecutive_count = 1
    print(indices_of_alignments)
    print(max_consecutive_count)
    print(index_longest_align)


def main():
    sequences = ["GTGACGTCACTCTGAGGATCCCCTGGGTGTGG",
                 "GTCAACTGCAACATGAGGAACATCGACAGGCCCAAGGTCTTCCT",
                 "GGATCCCCTGTCCTCTCTGTCACATA"]

    hash_table = construct_hash_table(sequences, 2)
    master_list = sequence_search("TGCAACAT", hash_table, 2)

    find_longest_match(master_list)

main()
