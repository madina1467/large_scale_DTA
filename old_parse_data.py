import csv

import pandas as pd
import numpy as np
import os
import json, pickle
from collections import OrderedDict
from rdkit import Chem
from rdkit.Chem import MolFromSmiles
import networkx as nx
from old_utils import *
from Bio import SeqIO

# def read_xml():

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


def create_dataset():
    fpath = '/Users/madina/Documents/research/denovo_moleculardesign/Custom-GraphDTA-with-BindingDB/large_scale_DTA/data/'
    fpath_save = '/Users/madina/Documents/research/denovo_moleculardesign/Custom-GraphDTA-with-BindingDB/large_scale_DTA/tt_batched/'

    ligands = csv.reader(open(fpath + "drugbank_db_smiles.csv"))
    fasta_sequences = SeqIO.parse(open(
        "data/uniprot-compressed_true_download_true_format_fasta_query__2A_20AND_2-2022.07.01-07.13.49.42.fasta"),
        'fasta')

    drugs = []
    prots = []

    header = next(ligands)
    for lg in ligands:
        lg = Chem.MolToSmiles(Chem.MolFromSmiles(lg[0]), isomericSmiles=True)
        drugs.append(lg)
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        prots.append(sequence)

    n_prots = len(prots) #20386
    n_drugs = len(drugs) #10920
    row_pairs = n_prots * n_drugs

    iterations = 48
    batch_size = row_pairs / iterations

    for iter in range(iterations):
        dataset = str(iter) + '_predict_' + str(batch_size)
        with open(fpath_save + dataset + '.csv', 'w') as f:
            f.write('Drug,Target,Y\n')  # Drug,Target,Y

            for prot_ind in range(n_prots):
                for drug_ind in range(n_drugs):
                    ls = []
                    ls += [drugs[drug_ind]]
                    ls += [prots[prot_ind]]
                    ls += [1]
                    f.write(','.join(map(str, ls)) + '\n')

        # print('\ndataset:', dataset)
        # print('train_fold:', len(train_fold))
        print('row_pairs:', row_pairs)
        # print('len(set(drugs)),len(set(prots)):', len(set(drugs)), len(set(prots)))
        # all_prots += list(set(prots))
        # compound_iso_smiles = set(drugs)






    # smile_graph = {}
    # for smile in drugs:
    #     g = smile_to_graph(smile)
    #     smile_graph[smile] = g
    #
    # processed_data_file = fpath + 'processed/' + dataset + '.pt'
    # if (not os.path.isfile(processed_data_file)):
    #     df = pd.read_csv(fpath + dataset + '.csv')
    #     test_drugs, test_prots,  test_Y = list(df['Drug']),list(df['Target']),list(df['Y'])
    #     XT = [seq_cat(t) for t in test_prots]
    #     test_drugs, test_prots = np.asarray(test_drugs), np.asarray(XT)
    #
    #     # make data PyTorch Geometric ready
    #     print('preparing ', dataset + '.pt in pytorch format!')
    #     # TODO !!!! y=[1]*row_pairs
    #     # test_data = TestbedDataset(root='data', dataset=dataset + '_test', xd=test_drugs, xt=test_prots,
    #     #                            y=[1] * row_pairs,
    #     #                            smile_graph=smile_graph)
    #     test_data = TestbedDataset(root=fpath, dataset=str(dataset + '.csv'), xd=test_drugs, xt=test_prots,
    #                                y=test_Y,
    #                                smile_graph=smile_graph)
    #     print(processed_data_file, ' and ', processed_data_file, ' have been created')
    # else:
    #     print(processed_data_file, ' and ', processed_data_file, ' are already created')


if __name__ == '__main__':
    create_dataset()
