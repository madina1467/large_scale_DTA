# from elasticsearch import Elasticsearch
# from elasticsearch import helpers
# import requests
#
# from rdkit.Chem import MolFromSmiles
#
#
# import os
#
# import numpy as np
#
# from sys import platform
# if platform == "darwin": # OS X
#     from large_scale_DTA.const import *
# else: # linux
#     from const import *
#
# import warnings
# import pandas as pd
# from pandas.core.common import SettingWithCopyWarning
# warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
#
# def addToDB(df, name, m_name):
#     es = Elasticsearch("http://localhost:9200")
#
#     col_names = list(df.columns)
#     actions = [
#             {
#                 col_names[0].lower(): row[col_names[0]],
#                 col_names[1].lower(): row[col_names[1]],
#                 col_names[3].lower().split("_", 2)[2]: row[col_names[3]],
#                 col_names[4].lower().split("_", 2)[2]: row[col_names[4]],
#                 col_names[5].lower().split("_", 2)[2]: row[col_names[5]],
#                 col_names[6].lower().split("_", 2)[2]: row[col_names[6]],
#                 col_names[7].lower().split("_", 2)[2]: row[col_names[7]],
#                 'model': m_name,
#             }
#             for index, row in df.iterrows()
#         ]
#
#     try:
#
#         response = helpers.bulk(es, actions, index='binding')
#         print("\nRESPONSE:", response)
#     except Exception as e:
#         print("\nERROR:", e)
#
# def prepare_processed_parquet():
#     print('Pre-processed data not found, doing pre-processing...')
#     df = pd.read_parquet(DATA_PATH + data_file + '.parquet')
#     test_drugs, test_prots,  test_Y = list(df['Drug']),list(df['Target']),list(df['Y'])
#     smile_graph = {}
#     for smile in set(test_drugs):
#         g = smile_to_graph(smile)
#         smile_graph[smile] = g
#
#     seq = dict()
#     XT = []
#     for t in test_prots:
#         if t not in seq.keys():
#             newseq = seq_cat(t)
#             XT.append(newseq)
#             seq[t] = newseq
#         else:
#             XT.append(seq[t])
#     test_drugs, xt = np.asarray(test_drugs), np.asarray(XT)
#     # test_data = TestbedDataset(root=DATA_PATH, dataset=str(data_file), xd=test_drugs, xt=xt, y=test_Y,
#     #                         smile_graph=smile_graph)
#
#     # test_loader = DataLoader(test_data, batch_size=TEST_BATCH_SIZE, shuffle=False)