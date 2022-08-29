from elasticsearch import Elasticsearch
from elasticsearch import helpers
import requests

from rdkit.Chem import MolFromSmiles
from large_scale_DTA.const import RESULTS_PATH, FINAL_RESULTS_PATH, TEST_RESULTS_PATH_TEST

from fastparquet import ParquetFile

import glob
import os

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

def createIndex():
    es = Elasticsearch("http://localhost:9200")

    if es.indices.exists(index="index"):
        body = "\"mappings\": {" \
        "\"_doc\": {" \
               "\"properties\": { \"Drug\": {" \
                "\"type\": \"text\"," \
                    "\"fields\": {" \
                        "\"keyword\": {" \
                            "\"ignore_above\": 256," \
                            "\"type\": \"keyword\" " \
                            "}" \
                    "}" \
               "}," \
               "\"GCNNet_kd\": {" \
                    "\"type\": \"float\"" \
               "}," \
               "\"Target\": {" \
                    "\"type\": \"text\"," \
                    "\"fields\": {" \
                        "\"keyword\": {" \
                            "\"ignore_above\": 256," \
                            "\"type\": \"keyword\"" \
                        "}" \
                    "}" \
               "}," \
               "\"GCNNet_ic50\": {" \
                    "\"type\": \"float\"" \
               "}," \
               "\"GATNet_ic50\": {" \
               "    \"type\": \"float\"" \
               "}," \
               "\"GAT_GCN_ic50\": {" \
               "    \"type\": \"float\"" \
               "}," \
               "\"GCNNet_ki\": {" \
               "    \"type\": \"float\"" \
               "}," \
               "\"GAT_GCN_ki\": {" \
               "    \"type\": \"float\"" \
                "}" \
            "}" \
        "}" \
    "}"
        es.indices.create(index="binding", body=body, params=None, headers=None)

def addToDB(df):
    es = Elasticsearch("http://localhost:9200")

    col_names = list(df.columns)
    # new_col_names = [column.replace("bdtdc_", "") for column in col_names]

    actions = []
    for index, row in df.iterrows():
        actions.append({
            col_names[i].replace("bdtdc_", ""): row[col_names[i]] for i in range(len(col_names)) if col_names[i].lower() != "y"
        })

        if index % 1 == 0: #13 times # 1356755
            try:
                response = helpers.bulk(es, actions, index='binding')
                print("\nRESPONSE:", response)
            except Exception as e:
                print("\nERROR:", e)
            actions = []


def delete_all():
    es = Elasticsearch("http://localhost:9200")
    es.indices.delete(index='binding', ignore=[400, 404])

    # def delete_action_gen():
    #     query_all = {
    #         'size': 10_000,
    #         'query': {
    #             'match_all': {}
    #         }
    #     }
    #     scanner = helpers.scan(es, index='binding',
    #                            query={'query': {'match_all': {}}, 'fields': []})
    #     for v in scanner:
    #         yield {'_op_type': 'delete',
    #                '_index': 'binding',
    #                '_type': v['_type'],
    #                '_id': v['_id'],
    #                }
    #
    # try:
    #     response = helpers.bulk(es, delete_action_gen())
    #     print("\nRESPONSE:", response)
    # except Exception as e:
    #     print("\nERROR:", e)


def appendToDB(df):
    es = Elasticsearch("http://localhost:9200")

    col_names = list(df.columns)
    # new_col_names = [column.replace("bdtdc_", "") for column in col_names]

    actions = []
    for index, row in df.iterrows():
        query = '{"query":{"bool":' \
                '{"must":[{"text":{"Drug":"' + row['Drug'] + '"}},' \
                         '{"text":{"Target":"' + row['Target'] + '"}}]' \
                '}},"from":0,"size":10,' \
                '"sort":[],"aggs":{}} '
        resp = es.search(index="binding", body=query)
        print("Got %d Hits:" % resp['hits']['total']['value'])
        results = resp['hits']['hits']
        for hit in resp['hits']['hits']:
            print(hit["_source"])

        # actions.append({
        #     new_col_names[i]: row[col_names[i]] for i in range(len(col_names)) if col_names[i].lower() != "y"
        # })
    #
    # try:
    #     response = helpers.bulk(es, actions, index='binding')
    #     print("\nRESPONSE:", response)
    # except Exception as e:
    #     print("\nERROR:", e)


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
        appendToDB(df)


def get_dataframe(file):
    filename = os.path.splitext(os.path.basename(file))[0]
    print('Transform file: ', filename)

    pf = ParquetFile(file)
    dataFrame = pf.to_pandas()
    return dataFrame
    # dataFrame.to_csv(FINAL_RESULTS_PATH + filename + '.csv', index=False)


if __name__ == '__main__':
    append_results_toDB()
    # delete_all()
    # add_results_toDB()

