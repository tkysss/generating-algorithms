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
import shutil
import argparse
import math

visualize_html_templates_dir = pathlib.Path(__file__).resolve().parent / 'visualize_html_templates'

visualization_root_dir = pathlib.Path('proof_visualization_output')

def get_graph_edges(graph):
    graph_edges = []
    for u in range(graph['n']):
        for v in range(u+1, graph['n']):
            edge_idx = (v*(v-1)//2+u)
            if edge_idx < len(graph['edges']) and graph['edges'][-1-edge_idx] == '1':
                graph_edges.append((u, v))
    return graph_edges

def get_name_for_graph_id(graph_id):
    return '{}_{}'.format(*graph_id)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('proof_export_dir_path')
    parser.add_argument('--visualization_dir_path')
    parser.add_argument('--overwrite_existing', action='store_true')

    args = parser.parse_args()

    proof_export_dir_path = pathlib.Path(args.proof_export_dir_path)
    if not args.visualization_dir_path:
        if args.proof_export_dir_path.endswith('.proof_export') and args.proof_export_dir_path.startswith('proof_export_output/'):
            p = args.proof_export_dir_path[len('proof_export_output/'):-len('.proof_export')] + '.proof_html_visualization'
            visualization_dir_path = visualization_root_dir / p
        else:
            visualization_dir_path = visualization_root_dir / (args.proof_export_dir_path + '.proof_html_visualization')
    else:
        visualization_dir_path = args.visualization_dir_path
    visualization_dir_path = pathlib.Path(visualization_dir_path)

    if visualization_dir_path.exists():
        if args.overwrite_existing:
            shutil.rmtree(visualization_dir_path.resolve())
        else:
            raise Exception('visualization_dir_path={!r} already exists, we do not want to overwrite anything'.format(visualization_dir_path))

    print('proof_export_dir_path={}'.format(proof_export_dir_path))
    print('visualization_dir_path={}'.format(visualization_dir_path))

    with open(proof_export_dir_path / 'proof.metadata.json', 'r') as f:
        proof_metadata = json.load(f)

    proof_visualizer = HTMLProofVisualizer(visualization_dir_path, proof_metadata)

    proof_visualizer.init()
    with gzip.open(proof_export_dir_path / 'proof.branching_rules.jsonl.gz', 'rt') as f:
        for line in f:
            d = json.loads(line)
            proof_visualizer.next(d)
            # break
    proof_visualizer.finalize()

class HTMLProofVisualizer:
    def __init__(self, visualization_dir_path, proof_metadata):
        self._visualization_dir_path = visualization_dir_path
        self._html_dir_path = visualization_dir_path / 'html'

        self._visualization_dir_path.mkdir(parents=True, exist_ok=True)
        self._html_dir_path.mkdir(parents=True, exist_ok=True)

        self._proof_metadata = proof_metadata
        self._graph_id_to_graph_num = proof_metadata['graph_id_to_graph_num']

        self._initial_graphs_index = []
        self._branching_rules_graphs_index = []
        self._all_graphs_index = []

        self._graph_id_to_parent_graph_id = {}
        self._graph_id_to_eliminated_parent_graph_id = defaultdict(set)

        self._algorithm_metadata = None

        self._jinja = jinja2.Environment(
            loader=jinja2.FileSystemLoader(visualize_html_templates_dir),
            autoescape=jinja2.select_autoescape()
        )

        self._next_cnt_so_far = 0

    def init(self):
        pass

    def finalize(self):
        self.finalize_algorithm_metadata()

        with open(self._html_dir_path / f'index.html', 'w') as f:
            template = self._jinja.get_template('index_html.j2')
            f.write(template.render({
                'proof_metadata': self._proof_metadata,
                'initial_graphs_index': self._initial_graphs_index,
                'branching_rules_graphs_index': self._branching_rules_graphs_index,
                'all_graphs_index': self._all_graphs_index,
                'algorithm_metadata': self._algorithm_metadata,
            }))

    def finalize_algorithm_metadata(self):
        graph_id_to_rule_tree_depth = {}
        def _compute_graph_id_to_rule_tree_depth(graph_id):
            if graph_id in graph_id_to_rule_tree_depth:
                return graph_id_to_rule_tree_depth[graph_id]
            d = None
            if self._graph_id_to_parent_graph_id[graph_id] == None:
                d = 0
            else:
                d = _compute_graph_id_to_rule_tree_depth(self._graph_id_to_parent_graph_id[graph_id]) + 1
            graph_id_to_rule_tree_depth[graph_id] = d
            return d
        for graph_id in self._graph_id_to_graph_num.keys():
            _compute_graph_id_to_rule_tree_depth(graph_id)
            assert graph_id in graph_id_to_rule_tree_depth

        graph_id_to_rule_walk_length = {}
        def _compute_graph_id_to_rule_walk_length(graph_id):
            if graph_id in graph_id_to_rule_walk_length:
                return graph_id_to_rule_walk_length[graph_id]
            d = None
            if self._graph_id_to_parent_graph_id[graph_id] == None:
                d = 0
            else:
                d = 0
                for parent_graph_id in self._graph_id_to_eliminated_parent_graph_id[graph_id]:
                    d = max(d, _compute_graph_id_to_rule_walk_length(parent_graph_id) + 1)
            graph_id_to_rule_walk_length[graph_id] = d
            return d
        for graph_id in self._graph_id_to_graph_num.keys():
            _compute_graph_id_to_rule_walk_length(graph_id)
            assert graph_id in graph_id_to_rule_walk_length

        rule_tree_depth_to_cnt = defaultdict(lambda: 0)
        for branching_rule in self._branching_rules_graphs_index:
            rule_tree_depth_to_cnt[graph_id_to_rule_tree_depth[branching_rule['graph_id']]] += 1

        rule_walk_length_to_cnt = defaultdict(lambda: 0)
        for branching_rule in self._branching_rules_graphs_index:
            rule_walk_length_to_cnt[graph_id_to_rule_walk_length[branching_rule['graph_id']]] += 1

        rule_tree_depth_rule_walk_length_to_cnt = defaultdict(lambda: 0)
        for branching_rule in self._branching_rules_graphs_index:
            rule_tree_depth_rule_walk_length_to_cnt[
                graph_id_to_rule_tree_depth[branching_rule['graph_id']] + graph_id_to_rule_walk_length[branching_rule['graph_id']]
            ] += 1

        self._algorithm_metadata['rule_tree_depth_to_cnt'] = rule_tree_depth_to_cnt
        self._algorithm_metadata['rule_walk_length_to_cnt'] = rule_walk_length_to_cnt
        self._algorithm_metadata['rule_tree_depth_rule_walk_length_to_cnt'] = rule_tree_depth_rule_walk_length_to_cnt

        def _d(d):
            if isinstance(d, dict):
                return {k: _d(v) for k, v in d.items()}
            return d

        self._algorithm_metadata['branching_rule_bf_rounded_hist_to_cnt'] = defaultdict(lambda: 0)
        for bf, cnt in self._algorithm_metadata['branching_rule_bf_to_cnt'].items():
            self._algorithm_metadata['branching_rule_bf_rounded_hist_to_cnt'][math.floor(bf*10)/10] += cnt

        self._algorithm_metadata['branching_rule_bf_rounded_hist_to_cnt'] = dict(self._algorithm_metadata['branching_rule_bf_rounded_hist_to_cnt'])

        self._algorithm_metadata['max_rule_tree_depth'] = max(self._algorithm_metadata['rule_tree_depth_to_cnt'].keys())
        self._algorithm_metadata['max_rule_walk_length'] = max(self._algorithm_metadata['rule_walk_length_to_cnt'].keys())
        self._algorithm_metadata['max_rule_tree_depth_rule_walk_length'] = max(self._algorithm_metadata['rule_tree_depth_rule_walk_length_to_cnt'].keys())
        self._algorithm_metadata['max_rule_size'] = max(self._algorithm_metadata['rule_size_to_cnt'].keys())
        self._algorithm_metadata['max_bf'] = max(self._algorithm_metadata['branching_rule_bf_to_cnt'].keys())

        self._algorithm_metadata['all_rules_cnt'] = self._algorithm_metadata['branching_rules_cnt'] + self._algorithm_metadata['reduction_rules_cnt']

        self._algorithm_metadata = _d(self._algorithm_metadata)

    def next(self, data):
        self._next_cnt_so_far += 1
        print('next progress {}/{}'.format(self._next_cnt_so_far, len(self._graph_id_to_graph_num)))
        _s = monotime()

        graph_id = data['graph']['id']
        branching_rule_data = {
            'proof_metadata': self._proof_metadata,
            'generation': data['generation'],
            'graph': {
                'graph_num': self._graph_id_to_graph_num[data['graph']['id']],
                'id': data['graph']['id'],
                'n': data['graph']['n'],
                'edges': get_graph_edges(data['graph']),
                'red_vertices': data['graph']['red_vs']
            },
            'parent_graph': {
                'graph_num': self._graph_id_to_graph_num[data['parent_graph']['id']],
                'id': data['parent_graph']['id'],
                'n': data['parent_graph']['n'],
                'edges': data['parent_graph']['edges'],
                'red_vertices': data['parent_graph']['red_vs'],
            } if data['parent_graph'] else None,
            'is_solved': data['is_solved'],
            'expansions': self.get_expansions_html_data(graph_id, data['expansions']),
            'rule': self.get_rule_html_data(data['rule']),
        }

        self.next_update_graph_index_data(branching_rule_data)
        self.next_update_algorithm_metadata(branching_rule_data)

        with open(self._html_dir_path / f'{graph_id}.html', 'w') as f:
            template = self._jinja.get_template('rule_html.j2')
            f.write(template.render(branching_rule_data))

    def next_update_graph_index_data(self, branching_rule_data):
        graph_id = branching_rule_data['graph']['id']
        graph_index_data = {
            'graph_id':  graph_id,
            'graph_num': self._graph_id_to_graph_num[graph_id],
            'rule_type': None
        }
        if branching_rule_data['rule']:
            graph_index_data['rule_type'] = branching_rule_data['rule']['rule_type']
            if graph_index_data['rule_type'] == 'branching_rule':
                graph_index_data['branching_rule_bf_rounded'] = branching_rule_data['rule']['bf_rounded']
            if graph_index_data['rule_type'] == 'reduction_rule':
                graph_index_data['reduction_rule_type'] = branching_rule_data['rule']['type']

        if branching_rule_data['generation'] == 1:
            self._initial_graphs_index.append(graph_index_data)
        if branching_rule_data['is_solved']:
            self._branching_rules_graphs_index.append(graph_index_data)
        self._all_graphs_index.append(graph_index_data)

    def next_update_algorithm_metadata(self, branching_rule_data):
        self._graph_id_to_parent_graph_id[branching_rule_data['graph']['id']] = branching_rule_data['parent_graph']['id'] if branching_rule_data['parent_graph'] else None
        for expansion in branching_rule_data['expansions']:
            assert expansion['expansion_result'] in ('eliminated', 'next', 'irrelevant')
            if expansion['expansion_result'] == 'eliminated':
                self._graph_id_to_eliminated_parent_graph_id[branching_rule_data['graph']['id']].add(expansion['result_data']['target_graph']['id'])

        if self._algorithm_metadata is None:
            self._algorithm_metadata = {
                'rule_size_to_cnt': defaultdict(lambda: 0),
                'branching_rules_cnt': 0,
                'branching_rule_type_to_cnt': defaultdict(lambda: 0),
                'branching_rule_bf_to_cnt': defaultdict(lambda: 0),
                'reduction_rules_cnt': 0,
                'reduction_rule_type_to_cnt': defaultdict(lambda: 0),
            }

        if branching_rule_data['rule']:
            self._algorithm_metadata['rule_size_to_cnt'][branching_rule_data['graph']['n']] += 1
            if branching_rule_data['rule']['rule_type'] == 'branching_rule':
                self._algorithm_metadata['branching_rule_bf_to_cnt'][branching_rule_data['rule']['bf_rounded']] += 1
                self._algorithm_metadata['branching_rules_cnt'] += 1
                subset_type_to_cnt = defaultdict(lambda: 0)
                for subset in branching_rule_data['rule']['subsets']:
                    subset_type_to_cnt[subset['type']] += 1
                assert subset_type_to_cnt['branch'] >= 1
                if subset_type_to_cnt['branch'] > 1:
                    self._algorithm_metadata['branching_rule_type_to_cnt']['branching'] += 1
                else:
                    self._algorithm_metadata['branching_rule_type_to_cnt']['reducing'] += 1
            elif branching_rule_data['rule']['rule_type'] == 'reduction_rule':
                self._algorithm_metadata['reduction_rule_type_to_cnt'][branching_rule_data['rule']['type']] += 1
                self._algorithm_metadata['reduction_rules_cnt'] += 1
            else:
                assert 0

    def get_expansions_html_data(self, graph_id, expansions):
        expansions_html_data = []
        for expansion_idx, expansion in enumerate(expansions, 1):
            expanding_vertex = expansion['graph']['n'] - 1
            expansion_graph_edges = get_graph_edges(expansion['graph'])

            expanding_graph_vertices = []
            for (u, v) in expansion_graph_edges:
                if expanding_vertex == u:
                    expanding_graph_vertices.append(v)
                if expanding_vertex == v:
                    expanding_graph_vertices.append(u)
            expanding_graph_vertices = sorted(set(expanding_graph_vertices))

            result = expansion['result']
            result_data = None
            if result == 'next' or result == 'eliminated':
                result_data = {
                    'target_graph': {
                        'graph_num': self._graph_id_to_graph_num[expansion['result_data']['target_graph']['id']],
                        'id': expansion['result_data']['target_graph']['id'],
                        'n': expansion['result_data']['target_graph']['n'],
                        'edges': get_graph_edges(expansion['result_data']['target_graph']),
                    },
                    'witnessing_isomorphism': expansion['result_data']['witnessing_isomorphism'],
                }

            expansion_html_data = {
                'expansion_graph': {
                    'graph_num': '{}_{}'.format(self._graph_id_to_graph_num[graph_id], expansion_idx),
                    'id': expansion['graph']['id'],
                    'n': expansion['graph']['n'],
                    'edges': get_graph_edges(expansion['graph']),
                },
                'expanding_graph_vertices': expanding_graph_vertices,
                'expansion_result': expansion['result'],
                'result_data': result_data,
            }
            expansions_html_data.append(expansion_html_data)

        return expansions_html_data

    def get_rule_html_data(self, rule):
        if not rule:
            return None

        if rule['type'] == 'handled_by_red_component_reduction':
            return {
                'rule_type': 'reduction_rule',
                'type': 'handled_by_red_component_reduction',
                'red_component1': rule['red_component1'],
                'red_component2': rule['red_component2'],
                'v': rule['v'],
            }
        elif rule['type'] == 'handled_by_red_star_reduction':
            return {
                'rule_type': 'reduction_rule',
                'type': 'handled_by_red_star_reduction',
                'kcenter_vertices': rule['kcenter_vertices'],
                'red_rays_vertices': rule['red_rays_vertices'],
            }
        elif rule['type'] == 'handled_by_vertex_cover_struction':
            return {
                'rule_type': 'reduction_rule',
                'type': 'handled_by_vertex_cover_struction',
                'v': rule['v'],
                'p': rule['p'],
                'q': rule['q'],
            }
        elif rule['type'] == 'branching_rule':
            return self.get_branching_rule_html_data(rule)
        else:
            assert 0, rule

    def get_branching_rule_html_data(self, rule):
        for subset in rule['subsets']:
            assert subset['type'] in ['branch', 'solution_but_adjusted', 'solution_but_dominated', 'solution_but_not_minimal', 'not_solution']
        return {
            'rule_type': 'branching_rule',
            'type': 'branching_rule',
            'subsets': rule['subsets'],
            'bf': rule['bf'],
            'bf_rounded': rule['bf_rounded'],
            'bv': rule['bv'],
        }

    # def get_branching_rule_html(self, rule):
    #     subset_rows_html = []
    #     for subset in rule['subsets']:
    #         # if subset['type'] == 'not_solution':
    #         #     subset_row_html = f"""
    #         #     <tr><td>{subset['subset']}</td><td>{subset['type']}</td><td>path: {subset['path']}</td></tr>
    #         #     """
    #         # elif subset['type'] == 'not_minimal_solution':
    #         #     subset_row_html = f"""
    #         #     <tr><td>{subset['subset']}</td><td>{subset['type']}</td><td>target: {subset['target']}</td></tr>
    #         #     """
    #         # elif subset['type'] == 'minimal_solution':
    #         #     subset_row_html = f"""
    #         #     <tr><td>{subset['subset']}</td><td>{subset['type']}</td><td></td></tr>
    #         #     """

    #         if subset['type'] == 'branch':
    #             subset_row_html = f"""
    #             <tr><td>{subset['subset']}</td><td>{subset['type']}</td><td></td></tr>
    #             """
    #         elif subset['type'] == 'solution_but_adjusted':
    #             subset_row_html = f"""
    #             <tr><td>{subset['subset']}</td><td>{subset['type']}</td><td>adjusted_by: {subset['adjusted_by']}</td></tr>
    #             """
    #         elif subset['type'] == 'solution_but_dominated':
    #             subset_row_html = f"""
    #             <tr><td>{subset['subset']}</td><td>{subset['type']}</td><td>dominated_by: {subset['dominated_by']}, red_vertices_subset: {subset['red_vertices_subset']}</td></tr>
    #             """
    #         elif subset['type'] == 'solution_but_not_minimal':
    #             subset_row_html = f"""
    #             <tr><td>{subset['subset']}</td><td>{subset['type']}</td><td>minimality_by: {subset['minimality_by']}</td></tr>
    #             """
    #         elif subset['type'] == 'not_solution':
    #             subset_row_html = f"""
    #             <tr><td>{subset['subset']}</td><td>{subset['type']}</td><td>path: {subset['path']}</td></tr>
    #             """
    #         else:
    #             assert 0, (subset, rule)

    #         subset_rows_html.append(subset_row_html)
    #     subset_rows_html = '\n'.join(subset_rows_html)

    #     rule_html = f"""
    #         <table>
    #         <tr><th>Subset</th><th>Type</th><th></th></tr>
    #         {subset_rows_html}
    #         </table>
    #     """

    #     return rule_html

    # def enqueue_generate_graph_image(self, name, n, edges, red_vertices):
    #     self._generate_graph_image_requests[name] = {
    #         'n': n,
    #         'edges': edges,
    #         'red_vertices': red_vertices,
    #     }

    # def generate_graph_images(self):
    #     with multiprocessing.Pool() as p:
    #         p.starmap(self.generate_graph_image, [
    #             (name, r['n'], r['edges'], r['red_vertices'])
    #             for name, r in self._generate_graph_image_requests.items()
    #         ])

    # def generate_graph_image(self, name, n, edges, red_vertices):
    #     _s = monotime()
    #     image_path = self._visualization_dir_path / f'images/{name}.png'
    #     if image_path.exists():
    #         return

    #     plt.clf()
    #     plt.axis('off')

    #     graph = nx.Graph()
    #     graph.add_nodes_from(range(n))

    #     node_color = []
    #     for v in range(n):
    #         node_color.append('red' if v in red_vertices else 'lightblue')

    #     graph.add_edges_from(edges)

    #     _ss = monotime()
    #     nx.draw_kamada_kawai(graph, with_labels=True, node_size=1000, width=2, font_size=20, node_color=node_color)

    #     plt.savefig(image_path)
    #     print('generate_graph_image took {:.4f}s'.format(monotime() - _s))

if __name__ == '__main__':
    main()
