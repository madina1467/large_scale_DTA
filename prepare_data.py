import os

import numpy as np
from const import *
import networkx as nx
from rdkit import Chem
import glob
from utils import TestbedDataset

import pandas as pd

seq_voc = "ABCDEFGHIKLMNOPQRSTUVWXYZ"
seq_dict = {v: (i + 1) for i, v in enumerate(seq_voc)}
seq_dict_len = len(seq_dict)
max_seq_len = 1000


def atom_features(atom):
    return np.array(one_of_k_encoding_unk(atom.GetSymbol(),
                                          ['C', 'N', 'O', 'S', 'F', 'Si', 'P', 'Cl', 'Br', 'Mg', 'Na', 'Ca', 'Fe', 'As',
                                           'Al', 'I', 'B', 'V', 'K', 'Tl', 'Yb', 'Sb', 'Sn', 'Ag', 'Pd', 'Co', 'Se',
                                           'Ti', 'Zn', 'H', 'Li', 'Ge', 'Cu', 'Au', 'Ni', 'Cd', 'In', 'Mn', 'Zr', 'Cr',
                                           'Pt', 'Hg', 'Pb', 'Unknown']) +
                    one_of_k_encoding(atom.GetDegree(), [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]) +
                    one_of_k_encoding_unk(atom.GetTotalNumHs(), [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]) +
                    one_of_k_encoding_unk(atom.GetImplicitValence(), [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]) +
                    [atom.GetIsAromatic()])


def one_of_k_encoding(x, allowable_set):
    if x not in allowable_set:
        raise Exception("input {0} not in allowable set{1}:".format(x, allowable_set))
    return list(map(lambda s: x == s, allowable_set))


def one_of_k_encoding_unk(x, allowable_set):
    """Maps inputs not in the allowable set to the last element."""
    if x not in allowable_set:
        x = allowable_set[-1]
    return list(map(lambda s: x == s, allowable_set))


def smile_to_graph(smile):
    mol = Chem.MolFromSmiles(smile)

    c_size = mol.GetNumAtoms()

    features = []
    for atom in mol.GetAtoms():
        feature = atom_features(atom)
        features.append(feature / sum(feature))

    edges = []
    for bond in mol.GetBonds():
        edges.append([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])
    g = nx.Graph(edges).to_directed()
    edge_index = []
    for e1, e2 in g.edges:
        edge_index.append([e1, e2])

    return c_size, features, edge_index


def seq_cat(prot):
    x = np.zeros(max_seq_len)
    for i, ch in enumerate(prot[:max_seq_len]):
        x[i] = seq_dict[ch]
    return x



if __name__ == '__main__':
    parquet_files = glob.glob(DATA_PATH_PARQUET + '*.parquet')

    for file in parquet_files:
        filename = os.path.splitext(os.path.basename(file))[0]
        if os.path.isfile(DATA_PATH_PARQUET_PROCESSED + filename + f'.pt'):
            print("\nPre-processed data => " + filename + ' -  data already created')

        else:
            print('\nPre-processed data not found, doing pre-processing...')
            df = pd.read_parquet(DATA_PATH_PARQUET + filename + '.parquet')
            test_drugs, test_prots, test_Y = list(df['Drug']), list(df['Target']), list(df['Y'])
            smile_graph = {}
            for smile in set(test_drugs):
                g = smile_to_graph(smile)
                smile_graph[smile] = g

            seq = dict()
            XT = []
            for t in test_prots:
                if t not in seq.keys():
                    newseq = seq_cat(t)
                    XT.append(newseq)
                    seq[t] = newseq
                else:
                    XT.append(seq[t])
            test_drugs, xt = np.asarray(test_drugs), np.asarray(XT)
            test_data = TestbedDataset(root=DATA_PATH_PARQUET, dataset=str(filename), xd=test_drugs, xt=xt, y=test_Y,
                                       smile_graph=smile_graph)
            print(filename + ".pt created in: ", DATA_PATH_PARQUET)