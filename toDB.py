from elasticsearch import Elasticsearch
from elasticsearch import helpers
import requests

from rdkit.Chem import MolFromSmiles
from large_scale_DTA.const import RESULTS_PATH, FINAL_RESULTS_PATH, TEST_RESULTS_PATH_TEST

from fastparquet import ParquetFile

import glob
import os
import base64

import numpy as np

from sys import platform

if platform == "darwin":  # OS X
    from large_scale_DTA.const import *
else:  # linux
    from const import *

import warnings
import pandas as pd
from pandas.core.common import SettingWithCopyWarning

warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)


def convertToBase64(message):
    # message_bytes = message.encode('ascii')
    # base64_bytes = base64.b64encode(message_bytes)
    # base64_message = base64_bytes.decode('ascii')
    # return base64.b64encode(message.encode('ascii'))
    # output = repr(base64.b64encode(bytes(message, "utf-8")))
    return base64.b64encode(bytes(message, "utf-8")).decode('ascii')


def createIndex():
    es = Elasticsearch("http://localhost:9200")

    che = {"settings": {"number_of_shards": 5, "number_of_replicas": 1},
           "mappings": {
               # "dynamic_templates": [
               #     {
               #         "strings": {
               #             "match_mapping_type": "string",
               #             "match": "*",
               #             "mapping": {
               #                 "type": "text",
               #                 "index": "not_analyzed"
               #             }
               #         }
               #     }
               # ],
               "properties": {
                   "Drug": {"type": "text"},
                   "Target": {"type": "text"},
                   "GCNNet_kd": {"type": "float"}, "GATNet_kd": {"type": "float"},
                   "GAT_GCN_kd": {"type": "float"}, "GCNNet_ic50": {"type": "float"},
                   "GATNet_ic50": {"type": "float"}, "GAT_GCN_ic50": {"type": "float"},
                   "GCNNet_ki": {"type": "float"}, "GATNet_ki": {"type": "float"},
                   "GAT_GCN_ki": {"type": "float"}, "MLTLE": {"type": "float"}
               }

           }
           }
    es.indices.create(index="binding", body=che)


def addToDB(df):
    es = Elasticsearch("http://localhost:9200")

    col_names = list(df.columns)
    # new_col_names = [column.replace("bdtdc_", "") for column in col_names]

    actions = []
    for index, row in df.iterrows():
        actions.append({
            col_names[i].replace("bdtdc_", ""): row[col_names[i]] for i in range(len(col_names)) if
            col_names[i].lower() != "y"
        })

        if index % 1 == 0:  # 13 times # 1356755
            try:
                response = helpers.bulk(es, actions, index='binding')
                print("\nRESPONSE:", response)
            except Exception as e:
                print("\nERROR:", e)
            actions = []


def delete_all():
    es = Elasticsearch("http://localhost:9200")
    es.indices.delete(index='binding', ignore=[400, 404])


def prepareFakeData():
    D1 = "CC(=O)NC(Cc1ccc2ccccc2c1)C(=O)NC(Cc1ccc(Cl)cc1)C(=O)NC(Cc1cccnc1)C(=O)NC(CO)C(=O)N(C)C(Cc1ccc(O)cc1)C(=O)NC(CC(N)=O)C(=O)NC(CC(C)C)C(=O)NC(CCCCNC(C)C)C(=O)N1CCCC1C(=O)NC(C)C(N)=O"
    T1 = "MFRRLTFAQLLFATVLGIAGGVYIFQPVFEQYAKDQKELKEKMQLVQESEEKKS"

    D2 = "CC(C)CC(NC(=O)CNC(=O)C(NC=O)C(C)C)C(=O)NC(C)C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(Cc1c[nH]c2ccccc12)C(=O)NC(CC(C)C)C(=O)NC(Cc1c[nH]c2ccccc12)C(=O)NC(CC(C)C)C(=O)NC(Cc1c[nH]c2ccccc12)C(=O)NC(CC(C)C)C(=O)NC(Cc1c[nH]c2ccccc12)C(=O)NCCO)C(C)C)C(C)C)C(C)C"
    T2 = "MVIMSEFSADPAGQGQGQQKPLRVGFYDIERTLGKGNFAVVKLARHRVTKTQVAIKIIDKTRLDSSNLEKIYREVQLMKLLNHPHIIKLYQVMETKDMLYIVTEFAKNGEMFDYLTSNGHLSENEARKKFWQILSAVEYCHDHHIVHRDLKTENLLLDGNMDIKLADFGFGNFYKSGEPLSTWCGSPPYAAPEVFEGKEYEGPQLDIWSLGVVLYVLVCGSLPFDGPNLPTLRQRVLEGRFRIPFFMSQDCESLIRRMLVVDPARRITIAQIRQHRWMRAEPCLPGPACPAFSAHSYTSNLGDYDEQALGIMQTLGVDRQRTVESLQNSSYNHFAAIYYLLLERLKEYRNAQCARPGPARQPRPRSSDLSGLEVPQEGLSTDPFRPALLCPQPQTLVQSVLQAEMDCELQSSLQWPLFFPVDASCSGVFRPRPVSPSSLLDTAISEEARQGPGLEEEQDTQESLPSSTGRRHTLAEVSTRLSPLTAPCIVVSPSTTASPAEGTSSDSCLTFSASKSPAGLSGTPATQGLLGACSPVRLASPFLGSQSATPVLQAQGGLGGAVLLPVSFQEGRRASDTSLTQGLKAFRQQLRKTTRTKGFLGLNKIKGLARQVCQVPASRASRGGLSPFHAPAQSPGLHGGAAGSREGWSLLEEVLEQQRLLQLQHHPAAAPGCSQAPQPAPAPFVIAPCDGPGAAPLPSTLLTSGLPLLPPPLLQTGASPVASAAQLLDTHLHIGTGPTALPAVPPPRLARLAPGCEPLGLLQGDCEMEDLMPCSLGTFVLVQ"

    newFakeDF = pd.DataFrame([[D1, T1, 3.5], [D1, T2, 3.5], [D2, T1, 3.5], [D2, T2, 3.5]],
                             columns=['Drug', 'Target', 'GCNNet_bdtdc_kd'])

    return newFakeDF


def appendToDB(df):
    es = Elasticsearch("http://localhost:9200")

    col_names = list(df.columns)
    # new_col_names = [column.replace("bdtdc_", "") for column in col_names]

    actions = []
    for index, row in df.iterrows():
        row['Drug'] = convertToBase64(row['Drug'])
        actions.append({
            col_names[i].replace("bdtdc_", ""): row[col_names[i]] for i in range(len(col_names)) if
            col_names[i].lower() != "y" and col_names[i].lower() != "test"
        })

    try:
        response = helpers.bulk(es, actions, index='binding')
        print("\nRESPONSE:", response)
        # print(actions)
    except Exception as e:
        print("\nERROR:", e)


def add_results_toDB():
    if not os.path.exists(TEST_RESULTS_PATH_TEST):
        os.makedirs(TEST_RESULTS_PATH_TEST)
        print("The new directory is created: ", TEST_RESULTS_PATH_TEST)

    parquet_files = glob.glob(TEST_RESULTS_PATH_TEST + '*_res_res.parquet')
    for file in parquet_files:
        df = get_dataframe(file)
        addToDB(df)


def append_results_toDB():
    if not os.path.exists(TEST_RESULTS_PATH_TEST):
        os.makedirs(TEST_RESULTS_PATH_TEST)
        print("The new directory is created: ", TEST_RESULTS_PATH_TEST)

    parquet_files = glob.glob(TEST_RESULTS_PATH_TEST + '*_res_res_fake_add.parquet')
    for file in parquet_files:
        df = get_dataframe(file)
        # appendToDB(df)
        # appendToDB()


def searchDF(df):
    es = Elasticsearch("http://localhost:9200")

    col_names = list(df.columns)
    # new_col_names = [column.replace("bdtdc_", "") for column in col_names]

    actions = []
    for index, row in df.iterrows():
        query = {
            "query": {
                "bool": {
                    "should": [
                        {
                            "query_string": {
                                "default_field": "Target",
                                "query": "*" + row['Target'] + "*"
                            }
                        },
                        {
                            "query_string": {
                                "default_field": "Drug",
                                "query": "*" + convertToBase64(row['Drug']) + "*"
                            }
                        }
                    ]
                }
            }
        }

        # print(query)

        # query_old = '{"query":{"bool":' \
        #         '{"must":[{"match":{"Drug":"' + row['Drug'] + '"}},' \
        #                   '{"match":{"Target":"' + row['Target'] + '"}}]' \
        #          '}},"from":0,"size":10, "sort":[],"aggs":{}} '
        resp = es.search(index="binding", body=query)
        print(resp)
        print("Got %d Hits:" % resp['hits']['total']['value'])
        results = resp['hits']['hits']
        for hit in resp['hits']['hits']:
            print(hit["_source"])


def search():
    if not os.path.exists(TEST_RESULTS_PATH_TEST):
        os.makedirs(TEST_RESULTS_PATH_TEST)
        print("The new directory is created: ", TEST_RESULTS_PATH_TEST)

    parquet_files = glob.glob(TEST_RESULTS_PATH_TEST + '*_res_res_fake_add.parquet')
    for file in parquet_files:
        df = get_dataframe(file)
        searchDF(df)


def get_dataframe(file):
    filename = os.path.splitext(os.path.basename(file))[0]
    print('Transform file: ', filename)

    pf = ParquetFile(file)
    dataFrame = pf.to_pandas()
    return dataFrame
    # dataFrame.to_csv(FINAL_RESULTS_PATH + filename + '.csv', index=False)


if __name__ == '__main__':
    # delete_all()
    # createIndex()
    fakeDF = prepareFakeData()
    # appendToDB(fakeDF)
    searchDF(fakeDF)

    # append_results_toDB()
    # search()



# {
#             "query": {
#                 "bool": {
#                     "should": [
#                         {
#                             "query_string": {
#                                 "default_field": "Target",
#                                 "query": "*" + row['Target'] + "*"
#                             }
#                         },
#                         {
#                             "query_string": {
#                                 "default_field": "Drug",
#                                 "query": "*" + convertToBase64(row['Drug']) + "*"
#                             }
#                         }
#                     ]
#                 }
#             }
#         }