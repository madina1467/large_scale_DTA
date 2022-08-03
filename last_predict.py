import os

import numpy as np
import torch
# from sys import platform
# if platform == "darwin": # OS X
#     from large_scale_DTA.const import *
# else: # linux
from const import *

from models.gcn import GCNNet
from models.gat_gcn import GAT_GCN

import networkx as nx
from rdkit import Chem

import glob
from tqdm.auto import tqdm

from utils import TestbedDataset
from torch_geometric.loader import DataLoader

import warnings
import pandas as pd
from pandas.core.common import SettingWithCopyWarning

warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)


def ged_models():
    md17 = 'BindingDB/GCNNet/model_GCNNet_bdtdc_kd.model'

    md31 = 'BindingDB/bd_tdc_latest_models/model_GCNNet_bdtdc_ki.model'  # 1
    md22 = 'BindingDB/GAT_GCN/model_GAT_GCN_bdtdc_ki.model'  # 2

    md24 = 'BindingDB/GCNNet/model_GCNNet_bindingdb_kd.model'  # 1
    md25 = 'BindingDB/GAT_GCN/model_GAT_GCN_bindingdb_kd.model'  # 2

    md26 = 'BindingDB/GCNNet/model_GCNNet_bindingdb_ki.model'  # 1
    md27 = 'BindingDB/GAT_GCN/model_GAT_GCN_bindingdb_ki.model'  # 2

    md28 = 'BindingDB/GAT_GCN/model_GAT_GCN_bindingdb_ic50.model'  # 2
    md29 = 'BindingDB/GCNNet/model_GCNNet_bindingdb_ic50.model'  # 1

    md30 = 'BindingDB/bd_tdc_latest_models/model_GCNNet_bdtdc_ic50.model'  # 1
    md23 = 'BindingDB/GAT_GCN/model_GAT_GCN_bdtdc_ic50.model'  # 2

    models_list = [md17, md22, md23, md25, md27, md28, md24, md26, md29, md30, md31]
    models_list = [MODELS_PATH + i for i in models_list]
    model_names = ['GCNNet_bdtdc_kd', 'GAT_GCN_bdtdc_ki', 'GAT_GCN_bdtdc_ic50', 'GAT_GCN_bindingdb_kd',
                   'GAT_GCN_bindingdb_ki', 'GAT_GCN_bindingdb_ic50',
                   'GCNNet_bindingdb_kd', 'GCNNet_bindingdb_ki', 'GCNNet_bindingdb_ic50', 'GCNNet_bdtdc_ic50',
                   'GCNNet_bdtdc_ki']
    model_arch = [GCNNet, GAT_GCN, GAT_GCN, GAT_GCN, GAT_GCN, GAT_GCN, GCNNet, GCNNet, GCNNet, GCNNet, GCNNet]
    m_name = ['GCNNet', 'GAT_GCN', 'GAT_GCN', 'GAT_GCN', 'GAT_GCN', 'GAT_GCN', 'GCNNet', 'GCNNet', 'GCNNet', 'GCNNet',
              'GCNNet']

    return models_list, model_names, model_arch, m_name


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


def seq_cat(prot):
    seq_voc = "ABCDEFGHIKLMNOPQRSTUVWXYZ"
    seq_dict = {v: (i + 1) for i, v in enumerate(seq_voc)}
    seq_dict_len = len(seq_dict)
    max_seq_len = 1000

    x = np.zeros(max_seq_len)
    for i, ch in enumerate(prot[:max_seq_len]):
        x[i] = seq_dict[ch]
    return x


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


def predicting(model, device, loader):
    model.eval()
    total_preds = torch.Tensor()
    # print('Make prediction for {} samples...'.format(len(loader.dataset)))
    with torch.no_grad():
        for data in tqdm(loader):
            data = data.to(device)
            output = model(data)
            total_preds = torch.cat((total_preds.to(device), output), 0)
    return total_preds.cpu().numpy().flatten()


def pred():
    models_list, model_names, model_arch, m_name = ged_models()

    if not os.path.exists(RESULTS_PATH):
        os.makedirs(RESULTS_PATH)
        print("The new directory is created: ", RESULTS_PATH)

    parquet_files = glob.glob(DATA_PATH_PARQUET + '*.parquet')
    print('Parquet data found, loading ...')
    # iter, iterations = 0, len(parquet_files)  # 48
    # batch_size = 4637815
    TEST_BATCH_SIZE = 10000

    # printProgressBar(iter, iterations, prefix='Progress:', suffix='Complete', length=iterations)
    # pbar = tqdm(total=batch_size, position=0, leave=True, ascii=True)


    for file in parquet_files:
        filename = os.path.splitext(os.path.basename(file))[0]
        if os.path.isfile(DATA_PATH_PARQUET_PROCESSED + filename + f'.pt'):
            print("Pre-processed data => " + filename + ' - found, loading ...')
            test_data = TestbedDataset(root=DATA_PATH_PARQUET, dataset=str(filename))
            test_loader = DataLoader(test_data, batch_size=TEST_BATCH_SIZE, shuffle=False)
        else:
            print("Sorry run prepare_data.py!!! Since no prepared data for " + filename + " file!!!")
            continue

            print('Pre-processed data not found, doing pre-processing...')
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
            test_loader = DataLoader(test_data, batch_size=TEST_BATCH_SIZE, shuffle=False)

        df_data = pd.read_parquet(DATA_PATH_PARQUET + filename + '.parquet')
        pred_df = df_data[df_data.index.isin(test_data.valid_id)]

        # pbar.update(1)

        for i in range(len(models_list)):
            device = torch.device('cuda:5' if torch.cuda.is_available() else "cpu")
            modeling = model_arch[i]
            model = modeling()
            model.to(device)
            model.load_state_dict(torch.load(models_list[i], map_location='cpu'))

            pred = predicting(model, device, test_loader)

            pred_df[model_names[i]] = pred

        pred_df.to_parquet(RESULTS_PATH + str(filename) + '_' + 'res.parquet', index=False)
        # printProgressBar(iter + 1, iterations, prefix='Progress:', suffix='Complete', length=iterations)
    # pbar.close()


if __name__ == '__main__':
    pred()
