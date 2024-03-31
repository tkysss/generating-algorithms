import sys
import json
import os
import networkx as nx
import matplotlib.pyplot as plt
import multiprocessing
import gzip
import itertools
from copy import deepcopy
from collections import defaultdict
import pathlib
from time import monotonic as monotime
import jinja2
import math

import argparse

def set_to_tuple(vs):
    return tuple(sorted(set(vs)))

def get_graph_edges(graph):
    graph_edges = []
    for u in range(graph['n']):
        for v in range(u+1, graph['n']):
            edge_idx = (v*(v-1)//2+u)
            if edge_idx < len(graph['edges']) and graph['edges'][-1-edge_idx] == '1':
                graph_edges.append((u, v))
    return graph_edges

def get_neighborhood_of_vertices(edges, vertices):
    vertices = set(vertices)
    N_vertices = set()

    for edge in edges:
        if len(set(edge) & vertices) == 1:
            N_vertices.add(edge[0])
            N_vertices.add(edge[1])

    return N_vertices - vertices

def tassert(text, t):
    tt = (text + '.' * 100)[:100]
    if t:
        pass
        # print(tt + 'OK')
    else:
        print(tt + 'FAILED')
        assert 0

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('proof_export_dir_path')

    args = parser.parse_args()

    proof_export_dir_path = pathlib.Path(args.proof_export_dir_path)

    with open(proof_export_dir_path / 'proof.metadata.json', 'r') as f:
        proof_metadata = json.load(f)

    proof_verifier = ProofVerifier(proof_metadata)

    proof_verifier.init()
    with gzip.open(proof_export_dir_path / 'proof.branching_rules.jsonl.gz', 'rt') as f:
        for idx, line in enumerate(f, 1):
            d = json.loads(line)
            proof_verifier.next(d, idx)
            # break
    proof_verifier.finalize()

class ProofVerifier:
    def __init__(self, proof_metadata):
        self._proof_metadata = proof_metadata
        self._graph_id_to_graph_num = proof_metadata['graph_id_to_graph_num']
        self._total_length = len(self._graph_id_to_graph_num)

    def init(self):
        pass

    def finalize(self):
        pass

    def next(self, data, idx):
        _s = monotime()

        self.verify_expansions(data)
        self.verify_red_vertices(data)
        if data['is_solved']:
            if data['rule']['type'] == 'branching_rule':
                self.verify_branching_rule(data)
            elif data['rule']['type'] in ('handled_by_red_component_reduction', 'handled_by_red_star_reduction', 'handled_by_vertex_cover_struction'):
                self.verify_reduction_rule(data)
            else:
                assert 0, data['rule']['type']

        print('{}/{}'.format(idx, self._total_length))

    def verify_expansions(self, data):
        graph = data['graph']

        # we create all connected supergraphs - all subsets of n vertices except the empty subset
        tassert('number of expansions', len(data['expansions']) == 2**data['graph']['n'] - 1)
        tassert('expansions all different', len(data['expansions']) == len(set([expansion['graph']['edges'] for expansion in data['expansions']])))

        for expansion in data['expansions']:

            expansion_graph = expansion['graph']

            tassert('expansion is supergraph', expansion_graph['edges'].endswith(graph['edges']) and expansion_graph['n'] == graph['n'] + 1)

            if expansion['result'] == 'next':
                target_graph = expansion['result_data']['target_graph']

                tassert('next target exists', target_graph['id'] in self._graph_id_to_graph_num)
                tassert('next target follows', self._graph_id_to_graph_num[graph['id']] < self._graph_id_to_graph_num[target_graph['id']])
            elif expansion['result'] == 'irrelevant':
                pass
            elif expansion['result'] == 'eliminated':
                target_graph = expansion['result_data']['target_graph']
                witnessing_isomorphism = expansion['result_data']['witnessing_isomorphism']

                tassert('eliminated target does not loop', self._graph_id_to_graph_num[graph['id']] > self._graph_id_to_graph_num[target_graph['id']])
                tassert('eliminated target isomorphism', self.check_isomorphism_subgraph(expansion_graph, target_graph, witnessing_isomorphism))
            else:
                assert 0, expansion

    def verify_red_vertices(self, data):
        red_vs = data['graph']['red_vs']
        if not red_vs:
            return

        red_vs_expansion_edges = set([(red_v, data['graph']['n']) for red_v in red_vs])
        for expansion in data['expansions']:
            if expansion['result'] != 'eliminated':
                expansion_edges = set(get_graph_edges(expansion['graph']))
                tassert('red vertices do not have edges in non eliminated expansions', len(expansion_edges & red_vs_expansion_edges) == 0)

    def verify_branching_rule(self, data):
        tassert('covering all non-empty subsets', len(set([tuple(sorted(branch['subset'])) for branch in data['rule']['subsets']])) == 2**data['graph']['n'] - 1)

        edges = set(get_graph_edges(data['graph']))

        subset_elimination_graph = {}
        subset_branches = set()
        for branch in data['rule']['subsets']:
            if branch['type'] == 'not_solution':
                continue
            elif branch['type'] == 'branch':
                subset_branches.add(set_to_tuple(branch['subset']))
            elif branch['type'] == 'solution_but_not_minimal':
                subset_elimination_graph[set_to_tuple(branch['subset'])] = set_to_tuple(branch['minimality_by'])
            elif branch['type'] == 'solution_but_adjusted':
                subset_elimination_graph[set_to_tuple(branch['subset'])] = set_to_tuple(branch['adjusted_by'])
            elif branch['type'] == 'solution_but_dominated':
                subset_elimination_graph[set_to_tuple(branch['subset'])] = set_to_tuple(branch['dominated_by'])
            else:
                assert 0, branch['type']

        subset_reaches_subset_branch = {}
        def _compute_subset_reaches_subset_branch(subset):
            if subset in subset_reaches_subset_branch:
                return subset_reaches_subset_branch[subset]
            if subset in subset_branches:
                r = True
            else:
                r = _compute_subset_reaches_subset_branch(subset_elimination_graph[subset])
            subset_reaches_subset_branch[subset] = r
            return r
        for branch in data['rule']['subsets']:
            if branch['type'] in ('solution_but_not_minimal', 'solution_but_adjusted', 'solution_but_dominated'):
                r = _compute_subset_reaches_subset_branch(set_to_tuple(branch['subset']))
                tassert('subset elimination order is correct', r == True)

        for branch in data['rule']['subsets']:
            if branch['type'] == 'branch':
                continue
            elif branch['type'] == 'not_solution':
                tassert('correct path len', len(branch['path']) == self._proof_metadata['dpvc_path_len'])
                tassert('path avoids solution', len(set(branch['path']) & set(branch['subset'])) == 0)
                path_edges = set()
                for i in range(len(branch['path'])-1):
                    u = branch['path'][i]
                    v = branch['path'][i+1]
                    if u > v:
                        u, v = v, u
                    path_edges.add((u, v))
                tassert('path exists in graph', len(path_edges & edges))
            elif branch['type'] == 'solution_but_not_minimal':
                tassert('check the minimality', set(branch['minimality_by']) <= set(branch['subset']))
            elif branch['type'] == 'solution_but_adjusted':
                tassert('check the adjustment', set(branch['adjusted_by']) <= set(branch['subset']))
            elif branch['type'] == 'solution_but_dominated':
                B = set(branch['subset'])
                R_star = set(branch['red_vertices_subset'])
                B_del = set(branch['subset']) - R_star
                B_d = set(branch['dominated_by'])

                N_H_prime_R_star = set()
                for edge in edges:
                    if len(set(edge) & B_del) != 0:
                        continue
                    if len(set(edge) & R_star) in (0, 2):
                        continue
                    N_H_prime_R_star.add(edge[0])
                    N_H_prime_R_star.add(edge[1])

                tassert('H[R*] is not bumpy', self.check_vertices_do_not_form_dpath_on_edges(R_star, edges))
                tassert("|R* \cap B| >= |N_H'(R*) \ R*|", len(R_star & B) >= len(N_H_prime_R_star - R_star))
                tassert("|N_H'(R*) \ R*| >= 1", len(N_H_prime_R_star - R_star) >= 1)
                tassert("B_d \subseteq (B \cup N_H'(R*)) \ R*", B_d <= ((B | N_H_prime_R_star) - R_star))
            else:
                assert 0, branch['type']

    def verify_reduction_rule(self, data):
        edges = set(get_graph_edges(data['graph']))
        red_vertices = set(data['graph']['red_vs'])

        rule = data['rule']

        if rule['type'] == 'handled_by_red_component_reduction':
            C_1 = rule['red_component1']
            C_2 = rule['red_component2']
            v = rule['v']

            tassert('C_1 is red', set(C_1) <= red_vertices)
            tassert('C_2 is red', set(C_2) <= red_vertices)
            tassert('C_1 is dpath free', self.check_vertices_do_not_form_dpath_on_edges(C_1, edges))
            tassert('C_2 is dpath free', self.check_vertices_do_not_form_dpath_on_edges(C_2, edges))

            N_C_1 = get_neighborhood_of_vertices(edges, C_1)
            N_C_2 = get_neighborhood_of_vertices(edges, C_2)
            tassert('N_C_1 is v', N_C_1 == {v})
            tassert('N_C_2 is v', N_C_1 == {v})
        elif rule['type'] == 'handled_by_red_star_reduction':
            C = rule['kcenter_vertices']
            L = rule['red_rays_vertices']
            d = self._proof_metadata['dpvc_path_len']

            tassert('|C| <= floor(d/2) - 1',  len(C) <= math.floor(d/2) - 1)
            tassert('|L| >= 2|C|', len(L) >= 2*len(C))
            tassert('L is red', set(L) <= red_vertices)
            for u in L:
                N_u = get_neighborhood_of_vertices(edges, [u])
                tassert('N_u = C', set(N_u) == set(C))
        elif rule['type'] == 'handled_by_vertex_cover_struction':
            v = rule['v']
            N_v = get_neighborhood_of_vertices(edges, [v])
            non_edges_in_N_v = len(N_v) * (len(N_v) - 1) / 2
            for edge in edges:
                if set(edge) <= N_v:
                    non_edges_in_N_v -= 1

            tassert('p is correct', rule['p'] == len(N_v))
            tassert('q is correct', rule['q'] == non_edges_in_N_v)
            tassert('p >= q', rule['p'] >= rule['q'])
        else:
            assert 0, rule['type']

    def check_isomorphism_subgraph(self, graph, subgraph, witnessing_isomorphism):
        # keys of the witnessing_isomorphism are the vertices of graph, values are the vertices of subgraph
        witnessing_isomorphism = {int(u): v for u, v in witnessing_isomorphism.items()}

        # first grab the edges of the graph
        graph_edges = get_graph_edges(graph)

        # grab only the induced edges - that is only edges whose both endpoints are being projected
        induced_graph_edges = [e for e in graph_edges if len(e & witnessing_isomorphism.keys()) == 2]

        # project the edges according to isomorphism
        projected_graph_edges = []
        for (u, v) in induced_graph_edges:
            fu = witnessing_isomorphism[u]
            fv = witnessing_isomorphism[v]
            if fu > fv:
                fu, fv = fv, fu
            projected_graph_edges.append((fu, fv))

        return sorted(projected_graph_edges) == sorted(get_graph_edges(subgraph))

    def check_vertices_do_not_form_dpath_on_edges(self, vertices, edges):
        if len(vertices) < self._proof_metadata['dpvc_path_len']:
            return True

        for path in itertools.permutations(sorted(vertices), self._proof_metadata['dpvc_path_len']):
            path_edges = set()
            for i in range(len(path)-1):
                u = path[i]
                v = path[i+1]
                if u > v:
                    u, v = v, u
                path_edges.add((u, v))
            if path_edges <= edges:
                return False
        return True

if __name__ == '__main__':
    main()
