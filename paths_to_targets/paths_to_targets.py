# -*- coding: utf-8 -*-
"""
Module finds all paths between a list of source body ids and targets body ids from Neuprint.
This is done by computing a breadth-first search starting with the the source body ids and
finding all downstream partners. This is repeated until the maximum search depth is
reached.

Examples:
    Execute paths_to_targests.py with a run file.  All parameters are set in tob.toml

        $ python paths_to_targests.py -f job.toml


    Execute paths_to_targests.py with a run file and override parameters command line
    argument.

        $ python paths_to_targests.py -f job.toml --roi_threshold 20 -p paths.csv

"""
from argparse import ArgumentParser
import ast
import re
from neuprint import Client
import numpy as np
import pandas as pd
import toml

class Graph:
    """
    Graph used to construct, store, and transform neuron connectivity graphs.

    Attributes:
        _client (Client): Neuprint API client.
        depth (int): Depth of search to build connectivity graph.
        graph (Dataframe): Rows are edges in connectivity graph.
        paths (Dataframe): Rows are paths in the format N0, W_1, N_1, ..., W_n, N_n.
        sources (list of int): List of source body IDs.
        targets (list of int): List of target body IDs.
    """


    def __init__(self, sources, targets, depth, verbose=False):
        """ Initialization method.
        Args:
            sources (list of int): List of source body IDs.
            targets (list of int): List of target body IDs.
            depth (int): Depth of search to build connectivity graph.

        Returns:
            None
        """
        self._client = None
        self.verbose = verbose
        self.depth = int(depth)
        self.graph = pd.DataFrame(columns=['SRC', 'DES', 'WEIGHT', 'LAYER'])
        self.paths = pd.DataFrame()
        self.sources = list(sources)
        self.targets = list(targets)


    def _verbose_print(self, msg):
        """ Verbocity print statement

        Args:
            msg (str): Message that should be printed for user
            verbose (bool):  If True print message

        Returns:
            None
        """
        if self.verbose:
            print(msg)


    def start_client(self, server, token):
        """ Creates a Neuprint API Client instance.

        Args:
            server (str): URL of Neuprint server.
            token (str): Authentication token for Neuprint.

        Returns:
            None
        """
        try:
            self._client = Client(server, token)
        except RuntimeError as error:
            raise RuntimeError("Issue instantiating Neuprint client:", error)


    def _query_downstream_partners(self, ids, threshold, hemibrain_neuron=True):
        """ Query Neuprint for downstream partners of a list of body IDs.

        Args:
            ids (list of int): List of body IDs to search from.
            threshold (int): Minimum connection strength between partners.
            hemibrain_neuron (bool): Specifies whether to search
                Hemibrain-Neuron of Hemibrain-Segment.

        Returns:
            Dataframe containing the connection pairs between query IDs
            and downstream partners.
        """
        query_template = (
            'WITH {IDS} AS SRC\n'
            'UNWIND SRC AS ID\n'
            'MATCH (n:`{HEMIBRAIN}`)-[w:ConnectsTo]->(m:`{HEMIBRAIN}`)\n'
            'WHERE n.bodyId = ID AND w.weight >= {THRESHOLD}\n'
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
        # print(query)
        results = self._client.fetch_custom(query)
        return results


    def _prune_graph(self, graph):
        """ Removes paths that do not include target bodies

        Starting from the top-most layer, remove all nodes from a layer that
        are not target bodies and do not connect to a body in the previous layer.

        Args:
            graph (Dataframe): Rows are edges of the connectivity graph.
            verbose (bool): When true print data regarding execution and
                progress.

        Returns:
            Dataframe with edges only edges that are part of pathways to targets
        """
        frontier = set(self.targets)
        pruned = pd.DataFrame(columns=['SRC', 'DES', 'WEIGHT', 'LAYER'])
        for layer in reversed(range(max(graph['LAYER']))):
            layer_edges = graph[graph['LAYER'] == layer + 1]
            lc1 = len(layer_edges)
            layer_edges = layer_edges[layer_edges.DES.isin(frontier)]
            pruned = pruned.append(layer_edges, ignore_index=True)
            frontier = set(layer_edges.SRC.values) | frontier
            self._verbose_print('Layer {LAYER}: {NUM_EDGES} edges'
                                .format(LAYER=layer + 1, NUM_EDGES=lc1))
            self._verbose_print('Layer {LAYER} Pruned: {NUM_EDGES} edges'
                                .format(LAYER=layer + 1, NUM_EDGES=len(layer_edges)))
        return pruned.reset_index(drop=True)


    def make_graph(self, threshold=1):
        """ Create connectivity graph of all downstream neurons

        Use breadth first search to find all bodies within N
        hops downstream of the source IDs.  Does not explore
        past nodes that have been visited or that are part of
        target set.

        Args:
            threshold (int, optional): Connection weight threshold between bodies
            verbose (bool, optional): Whether to print data

        Returns:
            Dataframe containing connectivity graph edges and connection weights.
        """
        frontier = set(self.sources)
        explored = set(self.targets)
        graph = pd.DataFrame(columns=['SRC', 'DES', 'WEIGHT', 'LAYER'])

        self._verbose_print('\nBuilding graph')
        for i in range(self.depth):
            self._verbose_print('Retrieving layer {LAYER} bodies'.format(LAYER=str(i+1)))

            layer_bodies = self._query_downstream_partners(frontier, threshold)
            layer_bodies['LAYER'] = i + 1
            explored = frontier | explored
            frontier = set(layer_bodies['DES'].values) - explored
            graph = graph.append(layer_bodies, ignore_index=True)

            self._verbose_print('Layer {LAYER}: {NUM} bodies'.format(NUM=len(frontier), LAYER=i+1))

        self._verbose_print('\nPruning graph of edges that do not lead to target bodies')
        self.graph = self._prune_graph(graph)
        return self.graph

    def graph_to_csv(self, file_name):
        """ Save graph to CSV

        Args:
            file_name (str): desired destination file for graph data

        Returns:
            None
        """
        self.graph.to_csv(file_name, index=False)


    def paths_to_csv(self, file_name):
        """ Save path to CSV

        Args:
            file_name (str): desired destination file for path data

        Returns:
            None
        """
        self.paths.to_csv(file_name, index=False)


    def compute_paths(self):
        """ Get all paths from from source bodies to targets.

        Find all permutations of edges that create pathways to target bodies.
        Starting from the second deepest layer, join layer N-1 DES column
        on layer N SRC column.

        Returns:

        """
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


def get_body_ids_from_roi(roi,
                          server,
                          token,
                          hemibrain_neuron=True,
                          pre_threshold=0,
                          post_threshold=0,
                          total_threshold=0):
    """

    Args:
        roi (str): Neuropil abreviation from Neuprint
        server (str): Neuprint server URL
        token (str): Neuprint auth token
        hemibrain_neuron (bool): Specifies whether to search
            Hemibrain-Neuron of Hemibrain-Segment.
        pre_threshold (int): Requires bodies meet threshold of presynaptic sites
            in the ROI
        post_threshold (int): Requires bodies meet threshold of postsynaptic sites
            in the ROI
        total_threshold (int): Requires bodies meet threshold of total synaptic sites
            in the ROI

    Returns:
        Dataframe containing body IDs that innervate the ROI

    """
    query_template = (
        'MATCH (n:`{HEMIBRAIN}`)\n'
        'WHERE n.{ROI}\n'
        'RETURN n.bodyId AS ID, n.name AS NAME, n.roiInfo AS ROIINFO'
    )

    # Start Neuprint python client
    client = Client(server, token)
    client.fetch_version()

    if hemibrain_neuron:
        query = query_template.format(HEMIBRAIN='hemibrain-Neuron', ROI=roi)
    else:
        query = query_template.format(HEMIBRAIN='hemibrain-Segment', ROI=roi)

    results = client.fetch_custom(query)
    results['ROIINFO'] = results['ROIINFO'].apply(ast.literal_eval)
    results['PRE'] = results['ROIINFO'].apply(lambda x: int(x[roi]['pre']))
    results['POST'] = results['ROIINFO'].apply(lambda x: int(x[roi]['post']))

    results = results[results['PRE'] + results['POST'] >= total_threshold]
    results = results[results['PRE'] >= pre_threshold]
    results = results[results['POST'] >= post_threshold]

    return results[['ID', 'NAME', 'PRE', 'POST']].reset_index()


def parse_arguments():
    """Creates dictionary of parameters from command line or from config file

    Returns:
        Dictionary of runtime parameters
    """
    args = dict()

    ##### Parse Commandline arguments #####
    parser = ArgumentParser()
    parser.add_argument('-f', '--file',
                        help='Input TOML containing execution parameters. Any parameters '
                             'that are manually specified '
                             'will override values in the file.')
    parser.add_argument('-S', '--server',
                        help='Neuprint server address')
    parser.add_argument('-T', '--token',
                        help='Neuprint API Token')
    parser.add_argument('-s', '--sources', nargs='+',
                        help='Comma separated list of body IDs where paths start.')
    parser.add_argument('-t', '--targets', nargs='+',
                        help='Comma separated list of body IDs for possible path destinations.')
    parser.add_argument('-r', '--roi_target',
                        help='Specify a Region Of Interest (ROI) from which to pull target bodies.')
    parser.add_argument('--roi_threshold', type=int,
                        help='Requires that bodies returned from the ROI have a combined total '
                             'presynaptic AND postsynaptic count greater than or equal to the '
                             'integer value specified.  Degault=0.')
    parser.add_argument('--roi_pre_threshold', type=int,
                        help='Requires that bodies returned from the ROI have a total '
                             'presynaptic count greater than or equal to the '
                             'integer value specified.  Degault=0.')
    parser.add_argument('--roi_post_threshold', type=int,
                        help='Requires that bodies returned from the ROI have a total '
                             'postsynaptic count greater than or equal to the '
                             'integer value specified.  Degault=0.')
    parser.add_argument('-d', '--depth', type=int,
                        help='Maximum length of paths.')
    parser.add_argument('-ct', '--connection_threshold', type=int,
                        help='Minimum connection strength between bodies in path.')
    parser.add_argument('-g', '--graph_file',
                        help='Output file location for the connection graph edges.')
    parser.add_argument('-p', '--path_file',
                        help='Output file location for paths')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Run script with runtime output.')
    cmd_args = parser.parse_args()

    ##### Read job file arguments #####
    if cmd_args.file:
        with open(cmd_args.file) as file:
            config = file.read()
            file_args = toml.loads(config)


    # Set run arguments to file values unless overridden by commandline arguemnts
    if cmd_args.server:
        args['server'] = cmd_args.server
    elif cmd_args.file and 'server' in file_args['client']:
        args['server'] = file_args['client']['server']
    else:
        raise RuntimeError('No server provided.')

    if cmd_args.token:
        args['token'] = cmd_args.token
    elif cmd_args.file and 'token' in file_args['client']:
        args['token'] = file_args['client']['token']
    else:
        raise RuntimeError('No token provided.')

    if cmd_args.sources:
        args['sources'] = cmd_args.sources
    elif cmd_args.file and 'sources' in file_args['job']:
        args['sources'] = file_args['job']['sources']
    else:
        raise RuntimeError('No source bodies provided.')

    sources = []
    if args['sources']:
        for arg in args['sources']:
            sources += re.findall(r'\d+', str(arg))
    args['sources'] = [int(source) for source in sources]

    if cmd_args.targets:
        args['targets'] = cmd_args.targets
    elif cmd_args.file and 'targets' in file_args['job']:
        args['targets'] = file_args['job']['targets']
    else:
        args['targets'] = None

    targets = []
    if args['targets']:
        for arg in args['targets']:
            targets += re.findall(r'\d+', str(arg))
    args['targets'] = [int(target) for target in targets]

    if cmd_args.roi_target:
        args['roi_target'] = cmd_args.roi_target
    elif cmd_args.file and 'roi_target' in file_args['job']:
        args['roi_target'] = file_args['job']['roi_target']
    else:
        args['roi_target'] = None

    # Test the either a Target or a Target ROI was created
    if not args['targets'] and not args['roi_target']:
        raise RuntimeError('No targets provided.')

    if cmd_args.roi_threshold:
        args['roi_threshold'] = cmd_args.roi_threshold
    elif cmd_args.file and 'roi_threshold' in file_args['job']:
        args['roi_threshold'] = file_args['job']['roi_threshold']
    else:
        args['roi_threshold'] = 0

    if cmd_args.roi_pre_threshold:
        args['roi_pre_threshold'] = cmd_args.roi_pre_threshold
    elif cmd_args.file and 'roi_pre_threshold' in file_args['job']:
        args['roi_pre_threshold'] = file_args['job']['roi_pre_threshold']
    else:
        args['roi_pre_threshold'] = 0

    if cmd_args.roi_post_threshold:
        args['roi_post_threshold'] = cmd_args.roi_post_threshold
    elif cmd_args.file and 'roi_post_threshold' in file_args['job']:
        args['roi_post_threshold'] = file_args['job']['roi_post_threshold']
    else:
        args['roi_post_threshold'] = 0

    if cmd_args.depth:
        args['depth'] = cmd_args.depth
    elif cmd_args.file and 'depth' in file_args['job']:
        args['depth'] = file_args['job']['depth']
    else:
        raise RuntimeError('No exploration depth provided')

    if cmd_args.connection_threshold:
        args['connection_threshold'] = cmd_args.connection_threshold
    elif cmd_args.file and 'connection_threshold' in file_args['job']:
        args['connection_threshold'] = file_args['job']['connection_threshold']
    else:
        raise RuntimeError('No connection threshold provided.')

    if cmd_args.graph_file:
        args['graph_file'] = cmd_args.graph_file
    elif cmd_args.file and 'graph_file' in file_args['job']:
        args['graph_file'] = file_args['job']['graph_file']
    else:
        args['graph_file'] = None

    if cmd_args.path_file:
        args['path_file'] = cmd_args.path_file
    elif cmd_args.file and 'path_file' in file_args['job']:
        args['path_file'] = file_args['job']['path_file']
    else:
        args['path_file'] = None

    if cmd_args.verbose:
        args['verbose'] = cmd_args.verbose
    elif cmd_args.file and 'path_file' in file_args['job']:
        args['verbose'] = file_args['job']['verbose']
    else:
        args['verbose'] = True

    return args


def main():
    """Main control flow if source_to_taget.py is called as a script.
    Returns:
        None
    """
    args = parse_arguments()
    sources = args['sources']
    targets = args['targets']
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
        raise RuntimeError('No targets found')

    if args['verbose']:
        print('Sources:', len(sources))
        print('Targets:', len(targets))

    graph = Graph(sources, targets, args['depth'], verbose=args['verbose'])
    graph.start_client(args['server'], args['token'])
    graph.make_graph(args['connection_threshold'])
    if args['graph_file']:
        graph.graph_to_csv(args['graph_file'])
    graph.compute_paths()
    if args['path_file']:
        graph.paths_to_csv(args['path_file'])


if __name__ == '__main__':
    main()
