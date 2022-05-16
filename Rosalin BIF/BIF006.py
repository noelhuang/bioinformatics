"""
Name: Noël Huang
Reg. nr: 1032641
Python 3 script to find the number of separate trees in a undirected graph.

Usage: CMD line: py BIF006.py <edge list format file>
Args:
    <edge list format file>: Edge list format in txt format. First line contains number
    of nodes and edges.
"""

import sys


def parse(file_name):
    """
    Parses an edge list format text file of a simple undirected graph.

    :param file_name: (str) name of the file to parse, must be an edge list format file. First row contains
    the number of vertices followed by number of edges.
    :return: data: list of list of edges

    Input example:
    2 3         # (vertices, edges)
    1 2         # (an edge)
    2 3         # (an edge)

    Output example: [[1, 2], [2, 3]]
    """
    with open(file_name) as f:
        data = f.read()
        data = data.split("\n")
        if data[-1] == "":
            del data[-1]
        for i in range(0, len(data)):
            data[i] = data[i].split(" ")
            for j in range(0, len(data[i])):
                data[i][j] = int(data[i][j])
    return data


def find_nodes(graph_data):
    """

    :param graph_data: List of an edge list format file. First element must contain number of vertices and edges.
    Example: [[2, 3], [1, 2], [2, 3]]

    :return: integer: number of nodes in the graph.
    """
    # Get number of nodes n and number of edges m
    nm = graph_data[0]
    number_of_nodes = nm[0]
    number_of_edges = nm[1]
    return number_of_nodes


def make_graph(graph_data):
    """
    Makes a cleaned graph from an edge list format.

    Args:
        graph_data: edge list format of a graph in text format. First line contains number of nodes and edges.

    Returns:
        Cleaned graph (number of nodes and edges is removed) in edge list format.

    """
    # Clean raw data by removing first line, which contains number of nodes and edges.
    del graph_data[0]
    return graph_data


def make_node_dict(number_of_nodes):
    """
    Makes a dictionary with all nodes in a graph as keys, with all of their values being False (bool).

    Args:
        number_of_nodes: (int) The number of nodes in a given graph.

    Returns: a dictionary with all nodes in a graph as keys, with all of their values being False (bool).

    The dictionary describes whether a node was visited already (True) or not (False).
    """
    # make dictionary of nodes (key) and a flag if visited (True) or unvisited (False)
    node_dic = {}
    for i in range(1, number_of_nodes + 1):
        node_dic[i] = False
    return node_dic


def explore(graph_data, some_node, node_dic):
    """
    Explores a graph from a start node (some_node) and finds all nodes that are accessible from this starting node.

    Args:
        graph_data: cleaned undirected graph as edge list format.
        some_node: A node in the graph, which is taken as the starting point for exploration of the graph.
        node_dic: Dictionary of all nodes in the graph as keys. The values indice whether this node was visited already
        or not (False/True).

    Returns: node_dic: Updated dictionary, which has value=True for all reachable nodes
     starting from the chosen start node (some_node).


    """
    a_node = some_node
    node_dic[a_node] = True # node_dic is a dictionary which contains flags for if a node was visited or not (True/False

    for edge in graph_data:
        if a_node in edge:
            index_a_node = edge.index(a_node)
            if index_a_node == 0:
                b_node = edge[1]
            elif index_a_node == 1:
                b_node = edge[0]
            if node_dic[b_node] == False:
                explore(graph_data, b_node, node_dic)
    return(node_dic)


def dfs(graph, node_dic):
    """
    Depth first search exploration of a graph to find all separate trees in the graph.

    Args:
        graph: cleaned graph (in edge list format)
        node_dic: Dictionary of all nodes in the graph as keys. The values indice whether this node was visited already
        or not (False/True).

    Returns: The number of separate trees in the given graph.

    """
    tree_count = 0
    for node in node_dic:
        if node_dic[node] == False:
            tree_count += 1
            explore(graph, node, node_dic)
    return(tree_count)


def main():
    file_name = sys.argv[1]
    graph_data = parse(file_name)
    number_of_nodes = find_nodes(graph_data)
    graph = make_graph(graph_data)
    node_dic = make_node_dict(number_of_nodes)
    print(dfs(graph, node_dic))


if __name__ == "__main__":
    main()

