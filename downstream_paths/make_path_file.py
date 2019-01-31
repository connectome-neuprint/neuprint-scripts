import numpy as np
import pandas as pd
import os
import re

nid_df = pd.read_csv('no_ids.csv')
SRC = [356818551,450034902,480029788]
DEST = list(nid_df[nid_df.pre >= 1].IDS.values)

count = 0

def get_downstream(id, layer, df, path, paths, cols):
    global count
    if id in DEST:
        path += [None] * (len(cols) - len(path))
        paths = paths.append(pd.Series(path, index=cols), ignore_index=True)
        count += 1
        return paths

    edges = df[df.SRC == id][df.Layer == layer]
    if len(edges.DES.values) != len(np.unique(edges.DES.values)):
        raise ValueError('There was a duplicate partner')
    ds_ids = list(edges.DES.values)
    i = 0
    for ds_id in ds_ids:
        i += 1
        if layer < 2:
            print(layer, "{n} of {total}".format(n=i, total=len(ds_ids)))
        if len(edges[edges.DES == ds_id]) > 1:
            raise ValueError('Duplicate ')
        paths = get_downstream(ds_id,
                               layer+1,
                               df,
                               path + [edges[edges.DES == ds_id].Weight.values[0],
                                       edges[edges.DES == ds_id].DES.values[0]],
                               paths,
                               cols)
    return paths


def main():
    print("Get all edges")
    df = pd.DataFrame(columns=['SRC','DES','Weight','Layer'])
    ls = os.listdir('fb')
    for file in [f for f in ls if 'edges_gte1_' in f]:
        dftmp = pd.read_csv(file, index_col=0)
        df = df.append(dftmp, ignore_index=True)

    # Prune Useless edges
    print("Prune useless edges")
    pruned = pd.DataFrame(columns=['SRC', 'DES', 'Weight', 'Layer'])
    prev = []
    for layer in reversed(range(max(df.Layer)+1)):
        dftmp = df[df.Layer == layer]
        dftmp = dftmp[dftmp.DES.isin(DEST+prev)]
        dftmp = dftmp[dftmp.SRC != dftmp.DES]
        pruned = pruned.append(dftmp, ignore_index=True)
        prev = list(np.unique(dftmp.SRC.values))
        print(layer, len(prev))
    pruned = pruned[pruned.Weight >= 5]

    prev = []
    for layer in range(max(df.Layer) + 1):
        dftmp = pruned[pruned.Layer == layer]
        dftmp = dftmp[dftmp.SRC.isin(SRC + prev)]
        dftmp = dftmp[dftmp.SRC != dftmp.DES]
        pruned = pruned[pruned.Layer != layer]
        pruned = pruned.append(dftmp, ignore_index=True)
        prev = list(np.unique(dftmp.DES.values))
        print(layer, len(prev))

    pruned = pruned.drop_duplicates()
    pruned.to_csv('fb/no_pruned_gte5.csv')

    # Create Graph File
    print("Making graph file")
    cols = ['SRC']
    for layer in range(max(df.Layer) + 1):
        cols += ['W{L}'.format(L=str(layer)), 'N{L}'.format(L=str(layer))]

    paths = pd.DataFrame(columns=cols)

    for bid in np.unique(pruned[pruned.Layer == 0].SRC):
        paths = get_downstream(bid,0, pruned, [bid], paths, cols)

    paths.to_csv('fb/paths_gte5_no.csv', index=False)



if __name__ == '__main__':
    main()
