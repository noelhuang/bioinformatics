"""
Name: Noël Huang
Reg. nr: 1032641
Python 3 script to find shortest path from a point to another in a simple directed graph.

Usage: CMD line: py BIF005.py <edge list format file> <number of starting vertex>
Args:
    <edge list format file>: Edge list format in txt format.
    <number of starting vertex>: Integer, which is the number of the starting vertex. Usually this is 1.
"""

import sys


def parse(file_name):
    """
    Parses an edge list format text file of a simple directed graph.

    :param file_name: (str) name of the file to parse, must be an edge list format file. First row contains
    the number of vertices followed by number of edges.
    :return: data: list of list of edges

    Input example:
    2 3         # (vertices, edges)
    1 2         # (an edge)
    2 3         # (an edge)

    Output example: [[1, 2], [2, 3]]
    """
    with open(file_name) as data:
        data = data.read()
        data = data.split("\n")
        if data[-1] == "":
            del data[-1]
        for i in range(0, len(data)):
            data[i] = data[i].split(" ")
            for j in range(0, len(data[i])):
                data[i][j] = int(data[i][j])
    # This data contains a list of lists with connected nodes (directed)
    return data


def bfs(data, starting_point):
    """
    Breadth-first search function for a simple directed graph. Outputs length of shortest path from starting point to i.

    :param data: A simple directed graph in edge list format, in .txt format
    :param starting_point: Number of the starting node (default 1).
    :return: String of integers separated by whitespace, representing the path length of starting_point to other nodes
    (sorted numerically from low to high).
    Negative integer (-1) indicates a node is not accessible from starting_point.
    """
    # Get information on total n, m. Where n = number of nodes, m= number of edges
    info = data.pop(0)
    n = info[0]
    m = info[1]
#    print(n, m)
#    data.sort()
    edge_list = data # Make list of edges for later use

    # Set distance all nodes to infinite
    node_dic = {}
    for i in range(1, n+1):
        node_dic[i] = "inf"
    # Set distance starting point to 0
    node_dic[starting_point] = 0
    queue = [starting_point] # This is our stack

    while len(queue) > 0:
        u = queue.pop(0) # u is a node number (starts with starting_point, but changes once stack is popped/pushed)
        for edge in edge_list:
            if edge[0] == u: # Starting point edge must be u
                if node_dic[edge[1]] == "inf": # Here we check if node was visited before, if not we add it to the queue
                    queue.append(edge[1])
                    node_dic[edge[1]] = node_dic[u] + 1 # Here we update the dictionary to add a unit distance to node u

    for key, value in node_dic.items():
        if value == "inf":
            node_dic[key] = -1
    values = [value for key, value in node_dic.items()]
    return values


def main():
    """
    Prints results, and writes them to output.txt.

    """
    file_name = sys.argv[1]
    number = int(sys.argv[2])
    data = parse(file_name)
    values = bfs(data, number)
    neat_values =" ".join(str(x) for x in values)
    print(neat_values)
    f = open("output.txt", "w")
    f.write(neat_values)


if __name__ == "__main__":
    main()

