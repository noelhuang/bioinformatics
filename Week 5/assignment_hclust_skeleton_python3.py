#!/usr/bin/env python3
"""
Author: Noël Huang
Student number:
Script to: 

Hints:
- write a function to obtain the Euclidean distance between two points.
- write a function to obtain the Distance matrix between all considered points
- write a function to obtain the correlation based distance between two points
- write a function to obtain the Distance matrix between all considered points
- write a function to implement h-clust
"""

# import statements go here 

from numpy import sqrt
import pprint


# function definitions

def csv_parser(lines):
    """Return list of point coordinates as [[x1,y1,z1,...],[x2,y2,z2,...]]
    
    lines: open file or list of lines. Expected format:
        The file has a single header line. 
        Each line contains the coordinates for one data point, starting
        with a label. A data point can be specified in arbitrary dimensions.

    Output: List of lists with the coordinates of the points.
    This function does not capture the labels of the data points. In case
    they are needed later, this function should be adjusted. 
    """

    data_points = []
    for line in lines:
        items = line.strip().split(",")
        try:  # will fail on header line in file
            data_points.append(list(map(float, items[1:])))  # skip label
        except ValueError:  # must be the header
            continue

    return data_points


class Clusterset(object):
    """Clusterset object describing a cluster

       # """

    def __init__(self, left=None, right=None, distance=0.0, ident=None):
        """ 
        ident: identifier of a leaf node (-1 for internal nodes)
        left: clusterset, left child of the node
        right: clusterset, right child of the node
        distance: distance between left and right children
        """
        self.left = left
        self.right = right
        self.ident = ident
        self.distance = distance


def print_tree(clust, labels=None, n=0):
    """Print a graphical representation of a clusterset object
    
    Input:
    clust: is a clusterset object. 
    n: integer to represent the deapth in the tree. n=0 is 
    the full tree. 
    labels: list of labels (strings) for elements in the tree
    """
    for i in range(n): print(' ', end='')  # end value required in python3
    if clust.ident < 0:  # negative ident means that this is branch
        print('-')
    else:  # positive ident means that this is an endpoint
        if labels == None:
            print(clust.ident)
        else:
            print(labels[clust.ident])
    # now print the right and left branches
    if clust.left != None:
        print_tree(clust.left, labels=labels, n=n + 1)
    if clust.right != None:
        print_tree(clust.right, labels=labels, n=n + 1)


def get_ordered_elements_distance(clust, list_of_elements=None, \
                                  list_of_dist=None):
    """Return ordered list of elements and distances from clusterset object
    
    Input:
    clust: Clusterset object
    list_of_elements: list with node identifiers. None is used to 
                       start the process.
    Output: 
    list_of_elements: list with node identifiers for the full tree
    list_of_dist: list of floats distance values 
    """
    if list_of_elements is None:
        list_of_elements = []
        list_of_dist = []
    if clust.ident < 0:  # negative ident means that internal
        list_of_elements.append('-')
        # append distance between the right and left branches
        list_of_dist.append(clust.distance)
    else:  # positive ident indicates leaf node (the dissim will be zero)
        list_of_elements.append(clust.ident)
        list_of_dist.append(clust.distance)
    if clust.left is not None:  # keep going left
        get_ordered_elements_distance(clust.left, \
                                      list_of_elements=list_of_elements, list_of_dist=list_of_dist)
    if clust.right is not None:  # keep going right
        get_ordered_elements_distance(clust.right, \
                                      list_of_elements=list_of_elements, list_of_dist=list_of_dist)

    return list_of_elements, list_of_dist


def cut_tree(list_of_elements, list_of_dist, h=1):
    """Return list of clusters from a list of elements 
    
    list_of_elements: list of node identifiers 
    list_of_dist: list of distance  values
    Both are obtained using get_ordered_elements_distance 
        on a clusterset object
    h: float with height to cut
    """
    clustlist = []
    cl = []
    # go through the elements and find brach point below the cut
    # add non internal nodes
    for i in range(len(list_of_elements)):
        if list_of_dist[i] is None:
            list_of_dist[i] = 0

        if list_of_dist[i] < h:
            if (list_of_elements[i] != '-'):
                cl.append(list_of_elements[i])

        if list_of_dist[i] >= h:
            if len(cl) > 0:
                clustlist.append(cl)
            cl = []
    clustlist.append(cl)
    return clustlist


def euclidean_distance(point_1, point_2):
    """
    Function to calculate the Euclidean distance between two points

    :param point_1: List of integers, representing the coordinates of point 1
    :param point_2: List of integers, representing the coordinates of point 2
    :return: Float, a scalar which represents the distance between point 1
    and point 2.
    """
    distance = sqrt(sum([(point_1[i] - point_2[i]) ** 2
                         for i in range(0, len(point_1))]))
    return distance


def euclidean_distance_matrix(data_points):
    number_of_points = len(data_points)
    matrix = [[None for i in range(0, number_of_points)]
              for i in range(0, number_of_points)]

    for i in range(0, number_of_points):
        for j in range(0, number_of_points):
            matrix[i][j] = euclidean_distance(data_points[i], data_points[j])

    return matrix


def correlation_distance(point_1, point_2):
    """
    Calculates Pearson correlation distance between two data points

    :param point_1: List of integers, representing the coordinates of point 1
    :param point_2: List of integers, representing the coordinates of point 2
    :return: Float, representing the pearson correlation distance between
    point 1 and point 2

    More info on Pearson Correlation Coefficient:
    https://www.wallstreetmojo.com/pearson-correlation-coefficient/
    https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
    """
    correlation = (len(point_1) * sum([point_1[i] * point_2[i]
                                       for i in range(0, len(point_1))])
                   - sum(point_1) * sum(point_2)) \
                  / (sqrt(len(point_1)
                          * sum([point_1[i] ** 2
                                 for i in range(0, len(point_1))])
                          - sum(point_1) ** 2)
                     * sqrt(len(point_2)
                            * sum([point_2[i] ** 2
                                   for i in range(0, len(point_2))])
                            - sum(point_2) ** 2))

    distance = 1 - correlation
    return distance


def correlation_distance_matrix(data_points):
    number_of_points = len(data_points)
    matrix = [[None for i in range(0, number_of_points)]
              for i in range(0, number_of_points)]

    for i in range(0, number_of_points):
        for j in range(0, number_of_points):
            matrix[i][j] = correlation_distance(data_points[i], data_points[j])

    return matrix


def distance_lookup(distance_dict, i, j):
    distance = distance_dict[i][j]
    return distance


def calc_new_cluster_distances(distance_dict, names_min_value):
    """
    Calculates the distance between a newly formed node and all other nodes.


    Uses the Pearson correlation distance to calculate the distance.
    """
    # Calculates the distance between a newly formed node, and the other data
    # points. Where a data point in this case
    # consists of a measurement of a gene over several timepoints.
    distance_dict_new = {}
    for i, rows in distance_dict.items():
        for j, columns in rows.items():
            distance_dict_new[j] = (distance_dict[names_min_value[0]][j] + \
                                    distance_dict[names_min_value[1]][j] - \
                                    distance_dict[names_min_value[0]][names_min_value[1]]) / 2

    return distance_dict_new


def hierarchical_clustering(distance_matrix):
    # Line 1 and 2, make dict with all nodes with index as keys
    node_dict = {}
    for index, row in enumerate(distance_matrix):
        node = Clusterset(left=None, right=None, distance=None, ident=index)
        node_dict[index] = node

    number_of_genes = len(node_dict)

    # Convert distance matrix into dict
    distance_dict = {}
    for row_index, row in enumerate(distance_matrix):
        distance_dict[row_index] = {}
        for column_index, column in enumerate(row):
            distance_dict[row_index][column_index] = distance_matrix[row_index][column_index]
    # print('dist dict', distance_dict)

    node_name_counter = number_of_genes
    # Line 3
    while len(node_dict) > 1:
        # for a in range(0, 1):
        # Line 4, find closest clusters
        min_value = 100
        coordinates_min_value = [None, None]
        for i, row in distance_dict.items():
            for j, column in row.items():
                if column < min_value:
                    if i != j:
                        min_value = column
                        names_min_value = [i, j]
        # print('min values', names_min_value)
        # for i, row in distance_dict.items():


        # Line 5, merge closest clusters
        left_node = node_dict[names_min_value[0]]

        right_node = node_dict[names_min_value[1]]

        new_node = Clusterset(left=left_node,
                              right=right_node,
                              ident=-1,
                              distance=distance_dict[names_min_value[0]][names_min_value[1]])  # TODO: this

        # Line 6, compute distance from new cluster to all other clusters
        new_distances = calc_new_cluster_distances(distance_dict, names_min_value)
        new_distances[node_name_counter] = 0
        distance_dict[node_name_counter] = new_distances # TODO: i add a new row here, but I need to make sure I also add the column at the end
        for row_name, row in distance_dict.items():
            if row_name != node_name_counter:
                row[node_name_counter] = new_distances[row_name]

        # Line 7, add new vertex C to the graph
        node_dict[node_name_counter] = new_node
        node_name_counter += 1
        node_dict.pop(names_min_value[0])
        node_dict.pop(names_min_value[1])

        # Line 8, remove rows and columns in the distance dict corresponding to
        # C1 and C2
        distance_dict.pop(names_min_value[0])
        distance_dict.pop(names_min_value[1])
        for row, dictio in distance_dict.items():
            dictio.pop(names_min_value[0])
            dictio.pop(names_min_value[1])
        print("Dict length:", len(node_dict))
    node_list = list(node_dict.values())
    return node_list[0]


def hierarchical_clustering_shortcut(distance_matrix):
    # Line 1 and 2, make dict with all nodes with index as keys
    node_dict = {}
    for index, row in enumerate(distance_matrix):
        node = Clusterset(left=None, right=None, distance=None, ident=index)
        node_dict[index] = node

    number_of_genes = len(node_dict)

    # Convert distance matrix into dict
    distance_dict = {}
    for row_index, row in enumerate(distance_matrix):
        distance_dict[row_index] = {}
        for column_index, column in enumerate(row):
            distance_dict[row_index][column_index] = distance_matrix[row_index][column_index]
    # print('dist dict', distance_dict)

    node_name_counter = number_of_genes
    # Line 3
    while len(node_dict) > 1:
        # for a in range(0, 1):
        # Line 4, find closest clusters
        min_value = 100
        coordinates_min_value = [None, None]
        for i, row in distance_dict.items():
            for j, column in row.items():
                if column < min_value:
                    if i != j:
                        min_value = column
                        names_min_value = [i, j]
                if min_value < 0.5:
                    break
            if min_value < 0.5:
                break
        # print('min values', names_min_value)
        # for i, row in distance_dict.items():


        # Line 5, merge closest clusters
        left_node = node_dict[names_min_value[0]]

        right_node = node_dict[names_min_value[1]]

        new_node = Clusterset(left=left_node,
                              right=right_node,
                              ident=-1,
                              distance=distance_dict[names_min_value[0]][names_min_value[1]])  

        # Line 6, compute distance from new cluster to all other clusters
        new_distances = calc_new_cluster_distances(distance_dict, names_min_value)
        new_distances[node_name_counter] = 0
        distance_dict[node_name_counter] = new_distances
        for row_name, row in distance_dict.items():
            if row_name != node_name_counter:
                row[node_name_counter] = new_distances[row_name]

        # Line 7, add new vertex C to the graph
        node_dict[node_name_counter] = new_node
        node_name_counter += 1
        node_dict.pop(names_min_value[0])
        node_dict.pop(names_min_value[1])

        # Line 8, remove rows and columns in the distance dict corresponding to
        # C1 and C2
        distance_dict.pop(names_min_value[0])
        distance_dict.pop(names_min_value[1])
        for row, dictio in distance_dict.items():
            dictio.pop(names_min_value[0])
            dictio.pop(names_min_value[1])
        print("Dict length:", len(node_dict))
    node_list = list(node_dict.values())
    return node_list[0]

def main():

    # Question 1
    print('Question 1:')
    with open("jp_fig10_1a.csv") as file:
        lines = file.readlines()
    parsed_data = csv_parser(lines)

    euclidean_matrix = euclidean_distance_matrix(parsed_data)
    for index, i in enumerate(euclidean_matrix):
        print(f"datapoint {index}", [round(number, 1) for number in i])

    # correlation_matrix = correlation_distance_matrix(parsed_data)
    # for index, i in enumerate(correlation_matrix):
    #     print(f"datapoint {index}", [round(number, 1) for number in i])

    master_node = hierarchical_clustering(euclidean_matrix)
    print_tree(master_node, labels=None, n=0)
    list_of_elements, list_of_dist = get_ordered_elements_distance(master_node)
    print('Question 1b:')
    print(cut_tree(list_of_elements, list_of_dist, h=3))

    print('Question 2a:')
    with open("example_gene_expression_four_genes.csv") as file:
        lines = file.readlines()
    parsed_data = csv_parser(lines)

    print('Euclidean matrix:')
    euclidean_matrix = euclidean_distance_matrix(parsed_data)
    for index, i in enumerate(euclidean_matrix):
        print(f"datapoint {index}", [round(number, 1) for number in i])

    print('\nQuestion 2b:\nCorrelation matrix')
    correlation_matrix = correlation_distance_matrix(parsed_data)
    for index, i in enumerate(correlation_matrix):
        print(f"datapoint {index}", [round(number, 1) for number in i])

    master_node = hierarchical_clustering(euclidean_matrix)
    print_tree(master_node, labels=None, n=0)
    list_of_elements, list_of_dist = get_ordered_elements_distance(master_node)
    print('Question 2c, euclidean clustering:')
    print(cut_tree(list_of_elements, list_of_dist, h=3))

    master_node = hierarchical_clustering(correlation_matrix)
    print_tree(master_node, labels=None, n=0)
    list_of_elements, list_of_dist = get_ordered_elements_distance(master_node)
    print('Question 2d, correlation clustering:')
    print(cut_tree(list_of_elements, list_of_dist, h=0.1))

    print('Question 7:')
    with open("proteomics_data.csv") as file:
        lines = file.readlines()
    parsed_data = csv_parser(lines)

    print('\nQuestion 7:\nCorrelation matrix')
    correlation_matrix = correlation_distance_matrix(parsed_data)
    # for index, i in enumerate(correlation_matrix):
    #     print(f"datapoint {index}", [round(number, 1) for number in i])

    ### HERE I AM USING THE SHORTCUT VERSION OF THE ALGORITHM ###
    master_node = hierarchical_clustering_shortcut(correlation_matrix)
    print_tree(master_node, labels=None, n=0)
    list_of_elements, list_of_dist = get_ordered_elements_distance(master_node)
    print('Question 7, correlation clustering with shortcut:')
    print(cut_tree(list_of_elements, list_of_dist, h=0.1))

if __name__ == "__main__":
    main()

