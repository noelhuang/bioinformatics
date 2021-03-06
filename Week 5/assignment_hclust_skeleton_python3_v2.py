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

    def __init__(self, left=None, right=None, distance=0.0, ident=None,
                 index=None, left_index=None, right_index=None):
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
        self.index = index
        self.left_index = left_index
        self.right_index = right_index


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


# def calc_new_cluster_distances(data_points, new_node):
#     """
#     Calculates the distance between a newly formed node and all other nodes.
#
#     :param data_points: list of lists of floats, can be seen as a 2D matrix
#     Where each row represents gene with measurements at different timepoints.
#     :param new_node: obj, A newly formed node, which is the combination of 2
#     existing nodes.
#     :return: List of floats, which represent the Pearson correlation distances
#     of the newly formed node to all the other nodes.
#
#     Uses the Pearson correlation distance to calculate the distance.
#     """
#     # Calculates the distance between a newly formed node, and the other data
#     # points. Where a data point in this case
#     # consists of a measurement of a gene over several timepoints.
#     distance_row = []
#     for data_point in data_points:
#         distance = (correlation_distance(data_point, new_node.left) +
#                     correlation_distance(data_point, new_node.right) -
#                     correlation_distance(new_node.left, new_node.right)) / 2
#         distance_row.append(distance)
#     return distance_row

def calc_new_distance_row(new_node, distance_matrix):
    distance_row = []
    for i in range(0, len(distance_matrix)):
        distance_c_star_c_1 = distance_matrix[i][new_node.left_index]
        distance_c_star_c_2 = distance_matrix[i][new_node.right_index]
        distance_c_1_c_2 = \
            distance_matrix[new_node.left_index][new_node.right_index]

        distance = (distance_c_star_c_1 + distance_c_star_c_2 -
                    distance_c_1_c_2) / 2
        distance_row.append(distance)
    distance_row.append(0)  # TODO: Think this is needed, check
    return distance_row


def hierarchical_clustering(distance_matrix):
    node_list = []
    for index, row in enumerate(distance_matrix):
        node = Clusterset(left=None, right=None, distance=None, ident=index,
                          index=index)
        node_list.append(node)

    while sum(isinstance(x, Clusterset) for x in node_list) > 2:
        # for a in range(0, 1):
        min_value = 100
        coordinates_min_value = [None, None]
        for j, column in enumerate(distance_matrix):
            for i, row in enumerate(column):
                if i != j:
                    if distance_matrix[i][j] < min_value:
                        min_value = distance_matrix[i][j]
                        coordinates_min_value = [i, j]

        # Find the correct nodes in the node list
        # for index, node in enumerate(node_list.copy()):
        #     if node_list[index] is not None:
        #         if node.index == coordinates_min_value[0]:
        #             index_sub_node_left = coordinates_min_value[0]
        #         if node.index == coordinates_min_value[1]:
        #             index_sub_node_right = coordinates_min_value[1]

        index_sub_node_left = coordinates_min_value[0]
        index_sub_node_right = coordinates_min_value[1]

        # Initiate new cluster with two children
        new_cluster = Clusterset(left=node_list[index_sub_node_left],
                                 right=node_list[index_sub_node_right],
                                 distance=distance_matrix[index_sub_node_left]
                                 [index_sub_node_right],
                                 ident=-1,
                                 index=len(distance_matrix),
                                 left_index=index_sub_node_left,
                                 right_index=index_sub_node_right)

        # Add new cluster to cluster list
        node_list.append(new_cluster)

        # Remove old nodes from cluster list
        node_list[index_sub_node_left] = None
        node_list[index_sub_node_right] = None

        # Set rows and columns corresponding to left and right subnodes to None
        # in the distance matrix
        distance_matrix[index_sub_node_left] = \
            [1000 for i in range(0, len(distance_matrix[index_sub_node_left]))]

        for i in range(len(distance_matrix)):
            distance_matrix[i][index_sub_node_right] = 1000

        # Add distance row and column for the newly formed cluster
        new_row = calc_new_distance_row(new_cluster, distance_matrix)
        distance_matrix.append(new_row)
        for i in range(0, len(distance_matrix)-1):
            distance_matrix[i].append(new_row[i])
        print(distance_matrix)

        print_tree(new_cluster, None, 0)


        last_node = node_list[-1]

        flag = True
        for x in new_row[0:-1]:
            if x < 3:
                flag = False

        if flag:
            node_list[-1] = None
        print("node list", node_list)
        print("newly formed node left index", new_cluster.left_index)
        print("newly formed node left clust obj", new_cluster.left)
        print("help this is count of clusters", sum(isinstance(x, Clusterset) for x in node_list))
    while None in node_list:
        node_list.remove(None)

    print("test", node_list)
    print(node_list[0].left_index)
    print(node_list[0].ident)
    print(node_list[0].index)
    print("last node", last_node.index)
    print("supertest")


def main():
    with open("jp_fig10_1a.csv") as file:
        lines = file.readlines()
    parsed_data = csv_parser(lines)

    euclidean_matrix = euclidean_distance_matrix(parsed_data)
    for index, i in enumerate(euclidean_matrix):
        print(f"datapoint {index + 1}", [round(number, 1) for number in i])

    correlation_matrix = correlation_distance_matrix(parsed_data)
    for index, i in enumerate(correlation_matrix):
        print(f"datapoint {index + 1}", [round(number, 1) for number in i])

    hierarchical_clustering(correlation_matrix)


if __name__ == "__main__":
    main()
    # the code below should produce the results necessary to answer
    # the questions. In other words, if we run your code, we should see 
    # the data that you used to answer the questions.
