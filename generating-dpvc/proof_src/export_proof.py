import sys
import json
import os
import networkx as nx
import matplotlib.pyplot as plt
import multiprocessing
import gzip
import itertools
from copy import deepcopy
from collections import defaultdict, deque
import pathlib
import scipy.optimize
import math
import shutil
import argparse

export_root_dir = pathlib.Path('proof_export_output')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('proof_input_file_path')
    parser.add_argument('--export_dir_path')
    parser.add_argument('--overwrite_existing', action='store_true')

    args = parser.parse_args()

    if not args.export_dir_path:
        if args.proof_input_file_path.endswith('.proof.out.gz') and args.proof_input_file_path.startswith('automated_improvement_output/'):
            p = pathlib.Path(args.proof_input_file_path[len('automated_improvement_output/'):-len('.proof.out.gz')] + '.proof_export')
            export_dir_path = export_root_dir / p
        else:
            export_dir_path = export_root_dir / (args.proof_input_file_path + '.proof_export')
    else:
        export_dir_path = args.export_dir_path
    export_dir_path = pathlib.Path(export_dir_path)

    if export_dir_path.exists():
        if args.overwrite_existing:
            shutil.rmtree(export_dir_path.resolve())
        else:
            raise Exception('export_dir_path={!r} already exists, we do not want to overwrite anything'.format(export_dir_path))

    export_dir_path.mkdir(parents=True, exist_ok=True)

    def proof():
        with gzip.open(args.proof_input_file_path, 'rt') as f:
            for line in f:
                d = json.loads(line)
                yield d
    proof_iter = iter(proof())

    ProofParser(PeekableProofIter(proof_iter), export_dir_path).parse_proof()

def get_graph_id(graph):
    return '{}_{}'.format(graph['n'], graph['edges'])

def get_canon_id(canon):
    return 'c{}_{}'.format(canon['c_n'], canon['c_edges'])

def get_graph_edges(graph):
    graph_edges = []
    for u in range(graph['n']):
        for v in range(u+1, graph['n']):
            edge_idx = (v*(v-1)//2+u)
            if edge_idx < len(graph['edges']) and graph['edges'][-1-edge_idx] == '1':
                graph_edges.append((u, v))
    return graph_edges

def get_graph_adj(graph):
    graph_edges = get_graph_edges(graph)
    graph_adj = defaultdict(list)
    for (u, v) in graph_edges:
        graph_adj[u].append(v)
        graph_adj[v].append(u)
    return dict(graph_adj)

def filter_out_graph_edges_with_vertex_v(graph_edges, v):
    return [
        e for e in graph_edges if v not in e
    ]

def get_graph_edges_mask_from_graph_edges(graph_edges):
    one_bits = set()
    for (u, v) in graph_edges:
        assert u < v
        edge_idx = (v*(v-1)//2+u)
        one_bits.add(edge_idx)
    max_edge_idx = max(one_bits)
    return ''.join(['1' if i in one_bits else '0' for i in range(max_edge_idx, -1, -1)])

def rename_vertex_u_to_v_in_graph_edges(graph_edges, u, v):
    new_graph_edges = []
    for (eu, ev) in graph_edges:
        if eu == u:
            eu = v
        if ev == u:
            ev = v
        if eu > ev:
            eu, ev = ev, eu
        assert eu < ev
        new_graph_edges.append((eu, ev))
    return new_graph_edges

def get_parent_graph_id(graph):
    pn = graph['n'] - 1
    pedges = graph['edges'][-(pn*(pn-1)//2):].lstrip('0') or '0'
    parent_graph_id = get_graph_id({
        'n': pn,
        'edges': pedges,
    })
    assert pedges != '0'
    return parent_graph_id

class PeekableProofIter:
    nopeek = object()

    def __init__(self,  proof_iter):
        self._proof_iter = proof_iter
        self._peek = self.nopeek


    def peek(self):
        if self._peek is self.nopeek:
            self._peek = next(self._proof_iter)
            self._proof_iter = itertools.chain([self._peek], self._proof_iter)
        return self._peek

    def next(self):
        self._peek = self.nopeek
        return next(self._proof_iter)

class ProofParser:
    def __init__(self, proof_iter, export_dir_path):
        self._proof_iter = proof_iter
        self._export_dir_path = export_dir_path
        self._proof_composer = None
        self._graph_num_counter = itertools.count(1)

    def peek(self):
        return self._proof_iter.peek()

    def next(self):
        r = self._proof_iter.next()
        # print(r)
        return r

    def parse_proof(self):
        self.parse_parameters()
        while True:
            try:
                self.peek()
            except StopIteration:
                break
            self.parse_generation()
        self._proof_composer.finalize()

    def parse_parameters(self):
        p = self.next()
        assert p['t'] == 'parameters', p

        dpvc_path_len = p['dpvc_path_len']
        dpvc_bf = p['dpvc_bf']
        generations = p['generations']

        self._proof_composer = ProofComposer(
            self._export_dir_path,
            dpvc_path_len,
            dpvc_bf,
            generations,
        )

    def parse_generation(self):
        p = self.next()
        assert p['t'] == 'generation_start', p

        generation = p['generation']
        print('parse_generation generation={}', generation)

        start_graphs = self.parse_current_bumpy_graphs()

        graph_id_to_filtered_expansions = self.parse_expansion()
        graph_id_to_rule = {}
        rule_graph_ids_in_order = []

        while True:
            p = self.peek()
            if p['t'] == 'generation_end':
                self.next()
                break
            round_graph_id_to_filtered_expansions, round_graph_id_to_rule, round_rule_graph_ids_in_order = self.parse_round()
            for graph_id, filtered_expansions in round_graph_id_to_filtered_expansions.items():
                graph_id_to_filtered_expansions[graph_id].extend(filtered_expansions)
            assert len(round_graph_id_to_rule.keys() & graph_id_to_rule.keys()) == 0
            graph_id_to_rule.update(round_graph_id_to_rule)
            rule_graph_ids_in_order.extend(round_rule_graph_ids_in_order)

        end_graphs, end_graphs_expansions = self.parse_current_bumpy_graphs_with_expansions()

        assert len(rule_graph_ids_in_order) == len(start_graphs) - len(end_graphs)
        graph_id_to_graph_num = {}
        for graph_id in rule_graph_ids_in_order:
            graph_id_to_graph_num[graph_id] = next(self._graph_num_counter)
        for graph in start_graphs:
            graph_id = get_graph_id(graph)
            if graph_id not in graph_id_to_graph_num:
                graph_id_to_graph_num[graph_id] = next(self._graph_num_counter)
        assert len(graph_id_to_graph_num) == len(start_graphs)
        for graph in start_graphs:
            graph['graph_num'] = graph_id_to_graph_num[get_graph_id(graph)]
        start_graphs.sort(key=lambda g: g['graph_num'])

        next_graphs = []
        try:
            if self.peek()['t'] == 'generation_next':
                self.next()
                next_graphs = self.parse_current_bumpy_graphs()
        except StopIteration:
            pass

        self._proof_composer.compose_generation(generation, start_graphs, graph_id_to_filtered_expansions, end_graphs, end_graphs_expansions, next_graphs, graph_id_to_rule)

    def parse_current_bumpy_graphs(self):
        current_bumpy_graphs = []

        p = self.next()
        assert p['t'] == 'current_bumpy_graphs_start'
        while True:
            p = self.next()
            assert p['t'] in ('current_bumpy_graph', 'current_bumpy_graphs_end')
            if p['t'] == 'current_bumpy_graphs_end':
                break

            current_bumpy_graphs.append(p['g'])

        return current_bumpy_graphs

    def parse_current_bumpy_graphs_with_expansions(self):
        current_bumpy_graphs = []
        current_bumpy_graphs_expansions = []

        p = self.next()
        assert p['t'] == 'current_bumpy_graphs_with_expansions_start'
        while True:
            p = self.next()
            if p['t'] == 'current_bumpy_graph':
                current_bumpy_graphs.append(p['g'])
                current_bumpy_graphs_expansions.append([])
            elif p['t'] == 'current_bumpy_graph_expansions_start':
                while True:
                    p = self.next()
                    assert p['t'] in ('current_bumpy_graph_expansion', 'current_bumpy_graph_expansions_end')
                    if p['t'] == 'current_bumpy_graph_expansions_end':
                        break
                    current_bumpy_graphs_expansions[-1].append(p['eg'])
            elif p['t'] == 'current_bumpy_graph_expansions_end':
                continue
            elif p['t'] == 'current_bumpy_graphs_with_expansions_end':
                break
            else:
                assert 0, p

        return current_bumpy_graphs, current_bumpy_graphs_expansions

    def parse_expansion(self):
        expand_gaffsf_match = []
        expand_gaffsf_cache = []

        p = self.next()
        assert p['t'] == 'expansion_start'
        while True:
            p = self.next()
            assert p['t'] in ('expansion_end', 'expand_gaffsf_match', 'expand_gaffsf_cache'), p
            if p['t'] == 'expansion_end':
                break
            if p['t'] == 'expand_gaffsf_match':
                expand_gaffsf_match.append(p)
            if p['t'] == 'expand_gaffsf_cache':
                expand_gaffsf_cache.append(p)

        return self.get_graph_id_to_filtered_expansions(expand_gaffsf_match, expand_gaffsf_cache)

    def parse_round(self):
        expand_gaffsf_match = []
        expand_gaffsf_cache = []

        rules = []

        p = self.next()
        assert p['t'] == 'round_start'
        while True:
            p = self.next()
            assert p['t'] in (
                'round_end',
                'expand_gaffsf_match', 'expand_gaffsf_cache',
                'handled_by_red_component_reduction', 'handled_by_red_star_reduction',
                'handled_by_vertex_cover_struction', 'handled_by_vertex_cover_dominance',
                'branching_rule',
                ), p
            if p['t'] == 'round_end':
                break
            if p['t'] == 'expand_gaffsf_match':
                expand_gaffsf_match.append(p)
            if p['t'] == 'expand_gaffsf_cache':
                expand_gaffsf_cache.append(p)
            if p['t'] in ('handled_by_red_component_reduction', 'handled_by_red_star_reduction', 'handled_by_vertex_cover_struction', 'handled_by_vertex_cover_dominance', 'branching_rule'):
                rules.append(p)

        return self.get_graph_id_to_filtered_expansions(expand_gaffsf_match, expand_gaffsf_cache), self.get_graph_id_to_rule(rules), [get_graph_id(p['g']) for p in rules]

    def get_graph_id_to_filtered_expansions(self, expand_gaffsf_match, expand_gaffsf_cache):
        graph_id_to_filtered_expansions = defaultdict(list)
        for p in expand_gaffsf_match:
            graph = p['g']
            expansion = p['eg']

            graph_id = get_graph_id(graph)
            graph_id_to_filtered_expansions[graph_id].append(p)

        for p in expand_gaffsf_cache:
            graph = p['g']
            expansion = p['eg']

            graph_id = get_graph_id(graph)

            graph_id_to_filtered_expansions[graph_id].append(p)

        return graph_id_to_filtered_expansions

    def get_graph_id_to_rule(self, rules):
        graph_id_to_rule = {}
        for p in rules:
            graph = p['g']
            graph_id = get_graph_id(graph)
            graph_id_to_rule[graph_id] = p
        return graph_id_to_rule

class ProofComposer:

    def __init__(self, export_dir_path, dpvc_path_len, dpvc_bf, generations):
        self._export_dir_path = export_dir_path
        self._proof_branching_rules_file = gzip.open(self._export_dir_path / 'proof.branching_rules.jsonl.gz', 'wt')
        self._proof_metadata_file = open(self._export_dir_path / 'proof.metadata.json', 'w')

        self._dpvc_path_len = dpvc_path_len
        self._dpvc_bf = dpvc_bf
        self._generations = generations

        self._solved_graph_id_to_graph = {}
        self._solved_canon_id_to_graph = {}
        self._graph_id_to_graph = {}
        self._canon_id_to_graph = {}

        self._graph_nums  = set()

        self._expansion_graph_id_to_match = {}

    def finalize(self):
        proof_metadata = {
            'dpvc_path_len': self._dpvc_path_len,
            'dpvc_bf': self._dpvc_bf,
            'graph_id_to_graph_num': {graph_id: g['graph_num'] for graph_id, g in self._graph_id_to_graph.items()},
            'canon_id_to_graph_num': {canon_id: g['graph_num'] for canon_id, g in self._canon_id_to_graph.items()},
        }
        self._proof_metadata_file.write(json.dumps(proof_metadata, indent=2))

    def compose_generation(self, generation, start_graphs, graph_id_to_filtered_expansions, end_graphs, end_graphs_expansions, next_graphs, graph_id_to_rule):
        for graph in start_graphs:
            graph_id = get_graph_id(graph)
            canon_id = get_canon_id(graph['canon'])
            assert graph_id not in self._graph_id_to_graph
            assert canon_id not in self._canon_id_to_graph
            self._graph_id_to_graph[graph_id] = graph
            self._canon_id_to_graph[canon_id] = graph
            assert graph['graph_num'] not in self._graph_nums
            self._graph_nums.add(graph['graph_num'])

        start_graph_ids = {get_graph_id(graph) for graph in start_graphs}
        end_graph_ids = {get_graph_id(graph) for graph in end_graphs}
        solved_graph_ids = start_graph_ids - end_graph_ids

        for graph_id in start_graph_ids:
            if graph_id in solved_graph_ids:
                assert graph_id in graph_id_to_rule
            else:
                assert graph_id not in graph_id_to_rule

        assert len(end_graphs) == len(end_graphs_expansions)
        end_graph_id_to_end_graph = {
            get_graph_id(end_graph): end_graph
            for end_graph in end_graphs
        }
        end_graph_id_to_graph_expansions = {
            get_graph_id(end_graph): end_graph_expansions
            for end_graph, end_graph_expansions in zip(end_graphs, end_graphs_expansions)
        }

        for solved_graph_id in solved_graph_ids:
            assert solved_graph_id not in self._solved_graph_id_to_graph
            graph = self._graph_id_to_graph[solved_graph_id]
            solved_canon_id = get_canon_id(graph['canon'])
            assert solved_canon_id not in self._solved_canon_id_to_graph
            self._solved_graph_id_to_graph[solved_graph_id] = graph
            self._solved_canon_id_to_graph[solved_canon_id] = graph

        for _idx, graph in enumerate(start_graphs, 1):
            print('composing graph {}/{}'.format(_idx, len(start_graphs)))
            graph_id = get_graph_id(graph)
            is_solved = graph_id in solved_graph_ids

            assert get_graph_edges_mask_from_graph_edges(get_graph_edges(graph)) == graph['edges']

            edges_count = graph['n'] * (graph['n'] - 1) // 2

            all_expansion_graphs = [{
                'n': graph['n'] + 1,
                'edges': bin(ei)[2:] + '0' * (edges_count - len(graph['edges'])) + graph['edges'],
            } for ei in range(1, 2**graph['n'])]

            all_expansion_graph_ids = [
                get_graph_id(eg) for eg in all_expansion_graphs
            ]

            expansion_graph_id_to_match = {}
            expansion_canon_id_to_match = {}
            expansion_graph_id_to_filtered_expansion = {}
            match_expansion_canon_id_to_expansion_graph_id = {}

            for filtered_expansion in graph_id_to_filtered_expansions[graph_id]:
                assert get_graph_id(filtered_expansion['g']) == graph_id
                expansion_graph_id = get_graph_id(filtered_expansion['eg'])
                expansion_graph_id_to_filtered_expansion[expansion_graph_id] = filtered_expansion
                assert expansion_graph_id in all_expansion_graph_ids


                assert filtered_expansion['t'] in ('expand_gaffsf_match', 'expand_gaffsf_cache')

                if filtered_expansion['t'] == 'expand_gaffsf_match':
                    expansion_canon_id = get_canon_id(filtered_expansion['eg']['canon'])
                    assert expansion_canon_id not in expansion_canon_id_to_match
                    match_canon = filtered_expansion['m_canon']
                    match_canon_id = get_canon_id(match_canon)
                    assert match_canon_id in self._solved_canon_id_to_graph
                    match_target_graph = self._solved_canon_id_to_graph[match_canon_id]

                    match_vertices = filtered_expansion['m_vs']

                    # We first project match_vertices to 0,1,2,...,len(m_vs)-1 and then project 0,1,2,...,len(m_vs)-1 to cv.
                    match_vertices_to_match_canon_isomorphism = {
                        mv: cv for mv, cv in zip(match_vertices, match_canon['c_relab'])
                    }

                    match_target_canon_to_match_target_vertices_isomorphism = {
                        cv: mv for mv, cv in enumerate(match_target_graph['canon']['c_relab'])
                    }

                    match_vertices_to_match_target_vertices_isomorphism = {
                        mv: match_target_canon_to_match_target_vertices_isomorphism[match_vertices_to_match_canon_isomorphism[mv]] for mv in match_vertices
                    }

                    match = {
                        'witnessing_isomorphism': match_vertices_to_match_target_vertices_isomorphism,
                        'target_graph': match_target_graph,
                    }
                    assert expansion_canon_id not in expansion_canon_id_to_match
                    expansion_canon_id_to_match[expansion_canon_id] = match
                    assert expansion_graph_id not in expansion_graph_id_to_match
                    assert expansion_graph_id not in self._expansion_graph_id_to_match
                    expansion_graph_id_to_match[expansion_graph_id] = match
                    self._expansion_graph_id_to_match[expansion_graph_id] = match
                    assert expansion_canon_id not in match_expansion_canon_id_to_expansion_graph_id
                    match_expansion_canon_id_to_expansion_graph_id[expansion_canon_id] = expansion_graph_id

                elif filtered_expansion['t'] == 'expand_gaffsf_cache':
                    expansion_graph_id = get_graph_id(filtered_expansion['eg'])
                    expansion_canon_id = get_canon_id(filtered_expansion['eg']['canon'])
                    assert expansion_canon_id in expansion_canon_id_to_match
                    match_filtered_expansion = expansion_graph_id_to_filtered_expansion[match_expansion_canon_id_to_expansion_graph_id[expansion_canon_id]]
                    assert expansion_canon_id == get_canon_id(match_filtered_expansion['eg']['canon'])

                    expansion_graph = filtered_expansion['eg']
                    match_expansion_graph = match_filtered_expansion['eg']

                    expansion_graph_canon_to_expansion_graph_vertices_isomorphism = {
                        cv: mv for mv, cv in enumerate(expansion_graph['canon']['c_relab'])
                    }

                    match_expansion_graph_vertices_to_expansion_graph_vertices_isomorphism = {
                        mv: expansion_graph_canon_to_expansion_graph_vertices_isomorphism[cv] for mv, cv in enumerate(match_expansion_graph['canon']['c_relab'])
                    }

                    match = expansion_canon_id_to_match[expansion_canon_id]

                    match = {
                        'witnessing_isomorphism': {
                            match_expansion_graph_vertices_to_expansion_graph_vertices_isomorphism[u]: v
                            for u, v in match['witnessing_isomorphism'].items()
                        },
                        'target_graph': match['target_graph'],
                    }
                    assert expansion_graph_id not in expansion_graph_id_to_match
                    assert expansion_graph_id not in self._expansion_graph_id_to_match
                    expansion_graph_id_to_match[expansion_graph_id] = match
                    self._expansion_graph_id_to_match[expansion_graph_id] = match
                else:
                    assert 0


            # Fill expansions of red vertices
            if graph['red_vs']:
                # If the graph has red_vertices, it cannot happen at the start of the first generation
                assert generation > 1
                # At the start of the generation, the last vertex added (from previous expansion) cannot be among the inherited red_vertices.
                assert graph['n'] - 1 not in graph['red_vs']
                parent_graph_id = get_parent_graph_id(graph)
                # The graph must have a parent graph we have already seen as this is not the first generation.
                assert parent_graph_id in self._graph_id_to_graph


                for eg in all_expansion_graphs:
                    expansion_graph_id = get_graph_id(eg)
                    eg_n = eg['n']
                    eg_edges = get_graph_edges(eg)
                    eg_adj = get_graph_adj(eg)

                    expands_red_vertices = len(set(graph['red_vs']) & set(eg_adj[eg_n - 1])) > 0

                    if not expands_red_vertices:
                        continue

                    assert expansion_graph_id not in expansion_graph_id_to_match
                    assert expansion_graph_id not in self._expansion_graph_id_to_match

                    # print('pre eg_edges', eg_edges)
                    sub_eg_edges = filter_out_graph_edges_with_vertex_v(eg_edges, eg_n - 2)
                    # print('post filter sub_eg_edges', sub_eg_edges)
                    sub_eg_edges = rename_vertex_u_to_v_in_graph_edges(sub_eg_edges, eg_n - 1, eg_n - 2)
                    # print('post rename sub_eg_edges', sub_eg_edges)
                    sub_eg_n = eg_n - 1

                    sub_eg = {
                        'n': sub_eg_n,
                        'edges': get_graph_edges_mask_from_graph_edges(sub_eg_edges)
                    }
                    sub_eg_id = get_graph_id(sub_eg)

                    # print('graph', graph)
                    # print('eg', eg)
                    # print('eg_n', eg_n)
                    # print('eg_edges', eg_edges)
                    # print('eg_adj', eg_adj)
                    # print('expands_red_vertices', expands_red_vertices)
                    # print('sub_eg_edges', sub_eg_edges)
                    # print('sub_eg_n', sub_eg_n)
                    # print('sub_eg', sub_eg)
                    # print('sub_eg_id', sub_eg_id)

                    assert sub_eg_id in self._expansion_graph_id_to_match

                    match = self._expansion_graph_id_to_match[sub_eg_id]
                    eg_witnessing_isomorphism = deepcopy(match['witnessing_isomorphism'])
                    eg_witnessing_isomorphism[eg_n - 1] = eg_witnessing_isomorphism[eg_n - 2]
                    eg_witnessing_isomorphism.pop(eg_n - 2)

                    match = {
                        'witnessing_isomorphism': eg_witnessing_isomorphism,
                        'target_graph': match['target_graph'],
                    }

                    expansion_graph_id_to_match[expansion_graph_id] = match
                    self._expansion_graph_id_to_match[expansion_graph_id] = match



            next_graph_canon_id_to_next_graph = {
                get_canon_id(next_graph['canon']): next_graph for next_graph in next_graphs
            }


            output_expansions = []

            # For each graph we always have to decide all its expansions.
            # Some of them are filtered in this generation (expansion_graph_id_to_match).
            # Others are filtered due to them expanding a red vertex (inherited red_vertices).
            # The remaining expansions were not eliminated and two situations occur:
            # Either the graph is solved, therefore the remaining expansions are irrelevant,
            # or the graph is not solved and each of the remaining expansions must be found in the next_graphs and will be solved in subsequent generations.
            for expansion_graph in all_expansion_graphs:
                expansion_graph_id = get_graph_id(expansion_graph)
                expansion_graph_result = None
                expansion_graph_result_data = None

                expansions_remaining_target_graph_id = None

                if is_solved:
                    if expansion_graph_id in expansion_graph_id_to_match:
                        target_graph = expansion_graph_id_to_match[expansion_graph_id]['target_graph']
                        expansion_graph_result = 'eliminated'
                        expansion_graph_result_data = {
                            'target_graph': {
                                'id': get_graph_id(target_graph),
                                'n': target_graph['n'],
                                'edges': target_graph['edges'],
                            },
                            'witnessing_isomorphism': expansion_graph_id_to_match[expansion_graph_id]['witnessing_isomorphism'],
                        }
                    else:
                        expansion_graph_result = 'irrelevant'
                else:
                    # If not solved, the graph must be among those at the end of generation.
                    assert graph_id in end_graph_id_to_graph_expansions
                    not_filtered_expansion_graph_id_to_expansion_graph = {
                        get_graph_id(eg): eg for eg in end_graph_id_to_graph_expansions[graph_id]
                    }
                    # The expansions must be either eliminated or not eliminated - not both.
                    assert (expansion_graph_id in expansion_graph_id_to_match) != (expansion_graph_id in not_filtered_expansion_graph_id_to_expansion_graph)
                    if expansion_graph_id in expansion_graph_id_to_match:
                        target_graph = expansion_graph_id_to_match[expansion_graph_id]['target_graph']
                        expansion_graph_result = 'eliminated'
                        expansion_graph_result_data = {
                            'target_graph': {
                                'id': get_graph_id(target_graph),
                                'n': target_graph['n'],
                                'edges': target_graph['edges'],
                            },
                            'witnessing_isomorphism': expansion_graph_id_to_match[expansion_graph_id]['witnessing_isomorphism'],
                        }
                    else:
                        e_expansion_graph = not_filtered_expansion_graph_id_to_expansion_graph[expansion_graph_id]
                        expansion_graph_canon_id = get_canon_id(e_expansion_graph['canon'])
                        assert expansion_graph_canon_id in next_graph_canon_id_to_next_graph
                        next_graph = next_graph_canon_id_to_next_graph[expansion_graph_canon_id]

                        next_graph_canon_to_next_graph_vertices_isomorphism = {
                            cv: mv for mv, cv in enumerate(next_graph['canon']['c_relab'])
                        }

                        expansion_graph_vertices_to_next_graph_vertices_isomorphism = {
                            mv: next_graph_canon_to_next_graph_vertices_isomorphism[cv] for mv, cv in enumerate(e_expansion_graph['canon']['c_relab'])
                        }

                        expansion_graph_result = 'next'
                        expansion_graph_result_data = {
                            'target_graph': {
                                'id': get_graph_id(next_graph),
                                'n': next_graph['n'],
                                'edges': next_graph['edges'],
                            },
                            'witnessing_isomorphism': expansion_graph_vertices_to_next_graph_vertices_isomorphism,
                        }

                assert expansion_graph_result is not None
                output_expansion = {
                    'graph': {
                        'id': expansion_graph_id,
                        'n': expansion_graph['n'],
                        'edges': expansion_graph['edges'],
                    },
                    'result': expansion_graph_result,
                    'result_data': expansion_graph_result_data,
                }
                output_expansions.append(output_expansion)


            for rv in graph['red_vs']:
                expanding_vertex = graph['n']
                for expansion in output_expansions:
                    expansion_adj = get_graph_adj(expansion['graph'])
                    if rv in expansion_adj[expanding_vertex]:
                        assert expansion['result'] == 'eliminated'



            parent_graph_id = None
            parent_graph = None
            if generation > 1:
                pn = graph['n'] - 1
                pedges = graph['edges'][-(pn*(pn-1)//2):].lstrip('0') or '0'
                parent_graph_id = get_graph_id({
                    'n': pn,
                    'edges': pedges,
                })

            if parent_graph_id is not None:
                assert parent_graph_id in self._graph_id_to_graph
                parent_graph = self._graph_id_to_graph[parent_graph_id]
                parent_graph['id'] = parent_graph_id

            is_solved = graph_id in solved_graph_ids

            # Update red_vertices to the end_graph.
            if not is_solved:
                graph['red_vs'] = end_graph_id_to_end_graph[graph_id]['red_vs']
                rule = None
            else:
                rule = graph_id_to_rule[graph_id]
                graph['red_vs'] = rule['g']['red_vs']
                rule = self.create_rule(graph, graph_id_to_rule[graph_id])

            graph['id'] = graph_id

            output_graph = {
                'generation': generation,
                'graph': graph,
                'parent_graph': parent_graph,
                'is_solved': is_solved,
                'rule': self.create_rule(graph, graph_id_to_rule[graph_id]) if is_solved else None,
                'expansions': output_expansions,
            }

            self._proof_branching_rules_file.write(json.dumps(output_graph) + '\n')

    def create_rule(self, graph, rule):
        # When we create a rule, first update the red_vertices.
        reduction_rules_t = ('handled_by_red_component_reduction', 'handled_by_red_star_reduction', 'handled_by_vertex_cover_struction', 'handled_by_vertex_cover_dominance')
        if rule['t'] == 'handled_by_red_component_reduction':
            return {
                'type': 'handled_by_red_component_reduction',
                'red_component1': rule['rc1'],
                'red_component2': rule['rc2'],
                'v': rule['v'],
            }
        elif rule['t'] == 'handled_by_red_star_reduction':
            return {
                'type': 'handled_by_red_star_reduction',
                'kcenter_vertices': rule['kcenter_vertices'],
                'red_rays_vertices': rule['red_rays_vertices'],
            }
        elif rule['t'] == 'handled_by_vertex_cover_struction':
            return {
                'type': 'handled_by_vertex_cover_struction',
                'v': rule['v'],
                'p': rule['p'],
                'q': rule['q'],
            }
        elif rule['t'] == 'handled_by_vertex_cover_dominance':
            assert 0, rule
        elif rule['t'] == 'branching_rule':
            return self.create_branching_rule(graph, rule)
        else:
            assert 0, rule

    def create_branching_rule(self, graph, rule):
        assert rule['t'] == 'branching_rule'
        subsets = []
        for s in range(1, 2**graph['n']):
            subset = [v for v in range(graph['n']) if ((1<<v)&s)]
            subsets.append(subset)

        graph_edges = set(get_graph_edges(graph))

        subset_to_data = {}

        minimal_solutions = rule['minimal_solutions']
        # The only minimal_solutions variant
        if not rule['dominance_free_solutions']:
            for subset in subsets:
                for minimal_solution in minimal_solutions:
                    if set(subset) > set(minimal_solution):
                        assert tuple(subset) not in subset_to_data
                        subset_to_data[tuple(subset)] = {
                            'type': 'solution_but_not_minimal',
                            'minimality_by': minimal_solution,
                        }
                        break

            for solution in minimal_solutions:
                subset_to_data[tuple(solution)] = {
                    'type': 'branch',
                }

            for subset in subsets:
                if tuple(subset) in subset_to_data:
                    continue
                path = self.get_path_in_subset_solution_for_graph(graph['n'], graph_edges, subset)
                assert path
                subset_to_data[tuple(subset)] = {
                    'type': 'not_solution',
                    'path': path,
                }
        else:
            dominance_free_solutions = rule['dominance_free_solutions']
            minimal_solution_to_dominating_solution = self.get_minimal_solution_to_dominating_solution(
                minimal_solutions,
                dominance_free_solutions,
                rule['solution_dominance_graph_adj'],
                rule['solution_dominance_graph_adj_red_vertices']
            )

            adjusted_solutions = rule['adjusted_solutions']
            dominance_free_solution_to_adjusting_solution = {}
            for solution in dominance_free_solutions:
                for a_solution in adjusted_solutions:
                    if set(a_solution) <= set(solution):
                        dominance_free_solution_to_adjusting_solution[tuple(solution)] = a_solution
                        break
            assert len(dominance_free_solutions) == len(dominance_free_solution_to_adjusting_solution)

            for solution in adjusted_solutions:
                solutions_being_adjusted = []
                for df_solution in dominance_free_solutions:
                    if set(solution) <= set(df_solution):
                        solutions_being_adjusted.append(df_solution)
                solutions_being_adjusted.sort()

                subset_to_data[tuple(solution)] = {
                    'type': 'branch',
                    # 'solutions_being_adjusted': solutions_being_adjusted,
                }

            for solution in dominance_free_solutions:
                a_solution = dominance_free_solution_to_adjusting_solution[tuple(solution)]
                if solution == a_solution:
                    subset_to_data[tuple(solution)] = {
                        'type': 'branch',
                    }
                else:
                    subset_to_data[tuple(solution)] = {
                        'type': 'solution_but_adjusted',
                        'adjusted_by': a_solution,
                    }

            for solution in minimal_solutions:
                if solution in dominance_free_solutions:
                    continue
                dominating_solution = minimal_solution_to_dominating_solution[tuple(solution)]
                subset_to_data[tuple(solution)] = {
                    'type': 'solution_but_dominated',
                    'dominated_by': dominating_solution['dominated_by'],
                    'red_vertices_subset': dominating_solution['red_vertices_subset'],
                }

            for subset in subsets:
                for minimal_solution in minimal_solutions:
                    if set(subset) > set(minimal_solution):
                        assert tuple(subset) not in subset_to_data
                        subset_to_data[tuple(subset)] = {
                            'type': 'solution_but_not_minimal',
                            'minimality_by': minimal_solution,
                        }
                        break

            for subset in subsets:
                if tuple(subset) not in subset_to_data:
                    path = self.get_path_in_subset_solution_for_graph(graph['n'], graph_edges, subset)
                    assert path
                    subset_to_data[tuple(subset)] = {
                        'type': 'not_solution',
                        'path': path,
                    }
                else:
                    assert subset_to_data[tuple(subset)]['type']

        for subset in subsets:
            assert tuple(subset) in subset_to_data

        output_subsets = []
        for subset in subsets:
            subset_data = subset_to_data[tuple(subset)]
            subset_data['subset'] = subset
            output_subsets.append(subset_data)

        branching_vector = []
        for subset in output_subsets:
            if subset['type'] == 'branch':
                branching_vector.append(len(subset['subset']))

        def branching_factor_f(x):
            return sum([(1/x)**b for b in branching_vector]) - 1

        root = scipy.optimize.root_scalar(branching_factor_f, x0=1, x1=self._dpvc_path_len, method='bisect', bracket=[1, self._dpvc_path_len])
        if root.converged:
            bf = root.root
            bf_rounded = math.ceil(bf * 10**4)/10**4
        else:
            bf = None
            bf_rounded = None

        return {
            'type': 'branching_rule',
            'subsets': output_subsets,
            'bf': bf,
            'bf_rounded': bf_rounded,
            'bv': sorted(branching_vector),
        }

    def get_path_in_subset_solution_for_graph(self, graph_n, graph_edges, subset):
        vertices = sorted(set(range(graph_n)) - set(subset))

        if len(vertices) < self._dpvc_path_len:
            return None

        for path_vertices in itertools.combinations(vertices, self._dpvc_path_len):
            for path_order in itertools.permutations(sorted(path_vertices)):
                is_path = True
                for u, v in zip(path_order[:-1], path_order[1:]):
                    if u > v:
                        u,v = v,u
                    assert u < v
                    if (u, v) not in graph_edges:
                        is_path = False
                        break
                if is_path:
                    return path_order
        return None

    def get_minimal_solution_to_dominating_solution(self, minimal_solutions, dominance_free_solutions, solution_dominance_graph_adj, solution_dominance_graph_adj_red_vertices):
        for solution in dominance_free_solutions:
            assert solution in minimal_solutions

        solution_tuple_to_solution_v = {
            tuple(solution): solution_v for solution_v, solution in enumerate(minimal_solutions)
        }


        closed = set()
        queue = deque()
        for solution in dominance_free_solutions:
            solution_v = solution_tuple_to_solution_v[tuple(solution)]
            closed.add(solution_v)
            queue.append(solution_v)

        solution_tuple_to_dominance = {}

        solution_dominance_graph_re_adj = defaultdict(list)
        solution_dominance_graph_re_adj_red_vertices = defaultdict(list)
        for u, (adj, adj_red_vertices) in enumerate(zip(solution_dominance_graph_adj, solution_dominance_graph_adj_red_vertices)):
            for v, v_red_vertices in zip(adj, adj_red_vertices):
                solution_dominance_graph_re_adj[v].append(u)
                solution_dominance_graph_re_adj_red_vertices[v].append(v_red_vertices)

        while len(queue):
            v = queue.popleft()

            for u, u_rv in zip(solution_dominance_graph_re_adj[v], solution_dominance_graph_re_adj_red_vertices[v]):
                if u not in closed:
                    closed.add(u)
                    queue.append(u)
                    solution_tuple_to_dominance[tuple(minimal_solutions[u])] = {
                        'dominated_by': minimal_solutions[v],
                        'red_vertices_subset': u_rv,
                    }

        for dominance_free_solution in dominance_free_solutions:
            assert tuple(dominance_free_solution) not in solution_tuple_to_dominance


        return solution_tuple_to_dominance

if __name__ == '__main__':
    main()
