__author__ ='Sorour E.Amiri'

import networkx as nx
import sys
import operator
import multiprocessing as mp


def generate_test_input():
    graph = {}
    graph['Source'] = {}
    graph['a'] = {}
    graph['b'] = {}
    graph['c'] = {}
    graph['d'] = {}
    graph['e'] = {}
    graph['f'] = {}
    graph['Target'] = {}
    graph['Source']['a'] = 3
    graph['Source']['d'] = 2
    graph['Source']['f'] = 10
    graph['a']['b'] = 2
    graph['a']['e'] = 7
    graph['d']['c'] = 4
    graph['f']['Target'] = 9
    graph['b']['c'] = 5
    graph['e']['Target'] = 8
    graph['c']['Target'] = 6
    return graph, 'Source', 'Target', 4


def make_dag(graph):
    dag = nx.DiGraph()
    for key1 in graph.keys():
        for key2 in graph[key1].keys():
            dag.add_weighted_edges_from([(key1, key2, graph[key1][key2])])
    return dag


def main(graph, source, target, max_length):
    alp_list = []
    dag = make_dag(graph)
    tree = nx.DiGraph()
    tree_nodes = {}
    tree_levels = {}
    level = 0
    predecessor = 'null'
    weight = 0
    c = 2
    tree_levels[level] = [source]
    tree_nodes[source] = {}
    tree_nodes[source][level] = (predecessor, weight)
    tree.add_node((source, level))

    for level in range(max_length + 1):
        tree_levels[level + 1] = []
        for nodes in tree_levels[level]:
            pi, p_weight = tree_nodes[nodes][level]
            for successor in dag.successors(nodes):
                # edge_weight = dag[nodes][successor]['weight']
                # tree.add_weighted_edges_from(((nodes, level), (successor, level + 1), edge_weight))
                # weight = p_weight + dag.edge[nodes][successor]['weight']
                weight = p_weight + dag.succ[nodes][successor]['weight']
                try:
                    curr_pi, curr_weight = tree_nodes[successor][level + 1]
                    if curr_weight < weight:
                        tree_nodes[successor][level + 1] = (nodes, weight)
                except KeyError:
                    try:
                        temp = tree_nodes[successor]
                    except KeyError:
                        tree_nodes[successor] = {}

                    tree_nodes[successor][level + 1] = (nodes, weight)
                    tree_levels[level + 1].append(successor)
        try:
            t_pi, t_weight = tree_nodes[target][level]
            alp_list.append(t_weight / float(level - c))
        except KeyError:
            alp_list.append(0)

    # alp_list = [lp_list[i] / float(i) for i in range(1, len(lp_list))]

    max_index, max_value = max(enumerate(alp_list), key=operator.itemgetter(1))

    # get the path
    # nx.all_simple_paths(tree, source=(source, 0), target=(target, max_index))
    pi, p_weight = tree_nodes[target][max_index]
    level = max_index - 1
    path = []
    while not level == 0:
        path.append(pi)
        pi, p_weight = tree_nodes[pi][level]
        level -= 1
    path.reverse()
    return path, max_value, max_index - c


if __name__ == '__main__':
    graph = sys.argv[1]
    source = sys.argv[2]
    target = sys.argv[3]
    max_length = sys.argv[4]
    # graph, source, target, max_length = generate_test_input()
    path, path_weight, path_length = main(graph, source, target, max_length)
    print(path, path_weight, path_length)