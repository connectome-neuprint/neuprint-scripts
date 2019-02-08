from argparse import ArgumentParser
import ast
import re
from neuprint import Client
import numpy as np
import pandas as pd
import toml


class Graph:
    def __init__(self, sources, targets, depth):
        self._client = None
        self.sources = list(sources)
        self.targets = list(targets)
        self.depth = int(depth)
        self.graph = pd.DataFrame(columns=['SRC', 'DES', 'WEIGHT', 'LAYER'])
        self.paths = pd.DataFrame()

    def start_client(self, server, token):
        self._client = Client(server, token)

    def _query_downstream_partners(self, ids, threshold, hemibrain_neuron=True):
        query_template = (
            'WITH {IDS} AS SRC\n'
            'MATCH (n:`{HEMIBRAIN}`)-[w:ConnectsTo]->(m:`{HEMIBRAIN}`)\n'
            'WHERE n.bodyId IN SRC AND w.weight >= {THRESHOLD}\n'
            'RETURN n.bodyId AS SRC, m.bodyId AS DES, w.weight AS WEIGHT'
        )

        if hemibrain_neuron:
            query = query_template.format(
                HEMIBRAIN='hemibrain-Neuron', IDS=list(ids), THRESHOLD=threshold
            )
        else:
            query = query_template.format(
                HEMIBRAIN='hemibrain-Segment', IDS=list(ids), THRESHOLD=threshold
            )
        results = self._client.fetch_custom(query)
        return results

    def _prune_graph(self, graph):
        frontier = set(self.targets)
        pruned = pd.DataFrame(columns=['SRC', 'DES', 'WEIGHT', 'LAYER'])
        for layer in reversed(range(max(graph['LAYER']))):
            layer_edges = graph[graph['LAYER'] == layer + 1]
            lc1 = len(layer_edges)
            layer_edges = layer_edges[layer_edges.DES.isin(frontier)]
            pruned = pruned.append(layer_edges, ignore_index=True)
            frontier = set(layer_edges.SRC.values) | frontier
            print(layer + 1, lc1, len(layer_edges))

        return pruned.reset_index()

    def make_graph(self, threshold=1):
        frontier = set(self.sources)
        explored = set(self.targets)
        graph = pd.DataFrame(columns=['SRC', 'DES', 'WEIGHT', 'LAYER'])

        for i in range(self.depth):
            layer_bodies = self._query_downstream_partners(frontier, threshold)
            layer_bodies['LAYER'] = i + 1
            explored = frontier | explored
            frontier = set(layer_bodies['DES'].values) - explored
            graph = graph.append(layer_bodies, ignore_index=True)

        self.graph = self._prune_graph(graph)
        return self.graph

    def graph_to_csv(self, file_name):
        self.graph.to_csv(file_name)

    def compute_paths(self):
        if self.graph.empty:
            raise Exception('No graph is not defined.')

        max_layer = self.graph['LAYER'].max()
        right = self.graph[self.graph['LAYER'] == max_layer][['SRC', 'WEIGHT', 'DES']]
        right = right.rename(
            index=str,
            columns={
                'SRC': 'N_{LAYER}'.format(LAYER=str(max_layer - 1)),
                'WEIGHT': 'W_{LAYER}'.format(LAYER=str(max_layer)),
                'DES': 'N_{LAYER}'.format(LAYER=str(max_layer)),
            },
        )

        for layer in reversed(range(max_layer - 1)):
            output_layer = 'N_{LAYER}'.format(LAYER=str(layer))
            connection_weight = 'W_{LAYER}'.format(LAYER=str(layer + 1))
            input_layer = 'N_{LAYER}'.format(LAYER=str(layer + 1))
            left = self.graph[self.graph['LAYER'] == layer + 1][
                ['SRC', 'WEIGHT', 'DES']
            ]
            left = left.rename(
                index=str,
                columns={
                    'SRC': output_layer,
                    'WEIGHT': connection_weight,
                    'DES': input_layer,
                },
            )
            right = pd.merge(left, right, on=input_layer, how='left')

        self.paths = right

    def paths_to_csv(self, file_name):
        self.paths.to_csv(file_name)


def get_body_ids_from_roi(roi,
                          server,
                          token,
                          hemibrain_neuron=True,
                          pre_threshold=0,
                          post_threshold=0,
                          total_threshold=0):
    client = Client(server, token)
    client.fetch_version()
    query_template = (
        'MATCH (n:`{HEMIBRAIN}`)\n'
        'WHERE n.{ROI}\n'
        'RETURN n.bodyId AS ID, n.name AS NAME, n.roiInfo AS ROIINFO'
    )

    if hemibrain_neuron:
        query = query_template.format(HEMIBRAIN='hemibrain-Neuron', ROI=roi)
    else:
        query = query_template.format(HEMIBRAIN='hemibrain-Segment', ROI=roi)

    results = client.fetch_custom(query)
    results['ROIINFO'] = results['ROIINFO'].apply(ast.literal_eval)
    results['PRE'] = results['ROIINFO'].apply(lambda x: int(x['FB']['pre']))
    results['POST'] = results['ROIINFO'].apply(lambda x: int(x['FB']['post']))

    results = results[results['PRE'] + results['POST'] >= total_threshold]
    print(pre_threshold)
    results = results[results['PRE'] >= pre_threshold]
    results = results[results['POST'] >= post_threshold]

    return results[['ID', 'NAME', 'PRE', 'POST']].reset_index()


def parse_arguments():
    args = dict()

    ##### Commandline arguments #####
    parser = ArgumentParser()
    parser.add_argument('-R', '--run_file')
    parser.add_argument('-S', '--server')
    parser.add_argument('-T', '--token')
    parser.add_argument('-s', '--sources', nargs='+')
    parser.add_argument('-t', '--targets', nargs='+')
    parser.add_argument('-r', '--roi_target')
    parser.add_argument('--roi_threshold', type=int)
    parser.add_argument('--roi_pre_threshold', type=int)
    parser.add_argument('--roi_post_threshold', type=int)
    parser.add_argument('-d', '--depth', type=int)
    parser.add_argument('-ct', '--connection_threshold', type=int)
    parser.add_argument('-g', '--graph_file')
    parser.add_argument('-p', '--path_file')
    cmd_args = parser.parse_args()

    ##### Run file arguments #####
    if cmd_args.run_file:
        with open(cmd_args.run_file) as file:
            config = file.read()
            file_args = toml.loads(config)

    if cmd_args.server:
        args['server'] = cmd_args.server
    elif cmd_args.run_file and 'server' in file_args['client']:
        args['server'] = file_args['client']['server']
    else:
        raise Exception('No server provided.')

    if cmd_args.token:
        args['token'] = cmd_args.token
    elif cmd_args.run_file and 'token' in file_args['client']:
        args['token'] = file_args['client']['token']
    else:
        raise Exception('No token provided.')

    if cmd_args.sources:
        args['sources'] = cmd_args.sources
    elif cmd_args.run_file and 'sources' in file_args['job']:
        args['sources'] = file_args['job']['sources']
    else:
        raise Exception('No source bodies provided.')

    if cmd_args.targets:
        args['targets'] = cmd_args.targets
    elif cmd_args.run_file and 'targets' in file_args['job']:
        args['targets'] = file_args['job']['targets']
    else:
        args['targets'] = None

    if cmd_args.roi_target:
        args['roi_target'] = cmd_args.roi_target
    elif cmd_args.run_file and 'roi_target' in file_args['job']:
        args['roi_target'] = file_args['job']['roi_target']
    else:
        args['roi_target'] = None

    # Test the either a Target or a Target ROI was created
    if not args['targets'] and not args['roi_target']:
        raise Exception('No targets provided.')

    if cmd_args.roi_threshold:
        args['roi_threshold'] = cmd_args.roi_threshold
    elif cmd_args.run_file and 'roi_threshold' in file_args['job']:
        args['roi_threshold'] = file_args['job']['roi_threshold']
    else:
        args['roi_threshold'] = 0

    if cmd_args.roi_pre_threshold:
        args['roi_pre_threshold'] = cmd_args.roi_pre_threshold
    elif cmd_args.run_file and 'roi_pre_threshold' in file_args['job']:
        args['roi_pre_threshold'] = file_args['job']['roi_pre_threshold']
    else:
        args['roi_pre_threshold'] = 0

    if cmd_args.roi_post_threshold:
        args['roi_post_threshold'] = cmd_args.roi_post_threshold
    elif cmd_args.run_file and 'roi_post_threshold' in file_args['job']:
        args['roi_post_threshold'] = file_args['job']['roi_post_threshold']
    else:
        args['roi_post_threshold'] = 0

    if cmd_args.depth:
        args['depth'] = cmd_args.depth
    elif cmd_args.run_file and 'depth' in file_args['job']:
        args['depth'] = file_args['job']['depth']
    else:
        raise Exception('No exploration depth provided')

    if cmd_args.connection_threshold:
        args['connection_threshold'] = cmd_args.connection_threshold
    elif cmd_args.run_file and 'connection_threshold' in file_args['job']:
        args['connection_threshold'] = file_args['job']['connection_threshold']
    else:
        raise Exception('No connection threshold provided.')

    if cmd_args.graph_file:
        args['graph_file'] = cmd_args.graph_file
    elif cmd_args.run_file and 'graph_file' in file_args['job']:
        args['graph_file'] = file_args['job']['graph_file']
    else:
        args['graph_file'] = None

    if cmd_args.path_file:
        args['path_file'] = cmd_args.path_file
    elif cmd_args.run_file and 'path_file' in file_args['job']:
        args['path_file'] = file_args['job']['path_file']
    else:
        args['path_file'] = None

    return args


def main():
    args = parse_arguments()

    # Parse Sources
    sources = []
    if args['sources']:
        for arg in args['sources']:
            sources += re.findall(r'\d+', str(arg))
    sources = [int(source) for source in sources]


    # Parse Targets
    targets = []
    if args['targets']:
        for arg in args['targets']:
            targets += re.findall(r'\d+', str(arg))
    targets = [int(target) for target in targets]
    if args['roi_target']:
        targets += list(
            get_body_ids_from_roi(
                args['roi_target'],
                args['server'],
                args['token'],
                pre_threshold=args['roi_pre_threshold'],
                post_threshold=args['roi_post_threshold'],
                total_threshold=args['roi_threshold'],
            )['ID'].values
        )
    # Remove dupilicates from target list
    targets = list(np.unique(targets))
    if not targets:
        raise Exception('No targets found')

    print('Sources:', len(sources))
    print('Targets:', len(targets))
    graph = Graph(sources, targets, args['depth'])
    graph.start_client(args['server'], args['token'])
    graph.make_graph(args['connection_threshold'])
    if args['graph_file']:
        graph.graph_to_csv(args['graph_file'])
    graph.compute_paths()
    if args['path_file']:
        graph.paths_to_csv(args['path_file'])


if __name__ == '__main__':
    main()
