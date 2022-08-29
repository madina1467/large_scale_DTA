import os
import pandas as pd
import tensorflow as tf
import numpy as np
import mltle as mlt
import shutil


def inference():
    CHUNK_SIZE = 50000
    MAX_BATCH_SIZE = 1024

    order = ['pKd', 'pKi', 'pIC50', 'is_active', 'qed']
    model = tf.keras.models.load_model('ResCNN1_pKd_pKi_pIC50.hdf5')



    filename = ''

    # path_data_predict = '0_drugbank_uniprot_4637815.csv'  # ~ 1923767
    path_data_predict = '/Users/madina/Documents/research/denovo_moleculardesign/Custom-GraphDTA-with-BindingDB/large_scale_DTA/data/batched/0_drugbank_uniprot_5.csv'  # ~ 1923767
    # path_result_save = '0_drugbank_uniprot_4637815_predict.csv'
    path_result_save = 'mlte_0_drugbank_uniprot_5.csv.csv'

    with open(path_result_save, 'w', encoding='utf-8') as f:
        f.write(f" ind, pKd , pKi, pIC50, isActive")

    # %%time

    DATA_CHUNKS = pd.read_csv(path_data_predict, chunksize=CHUNK_SIZE, usecols=['Drug', 'Target'])

    for data in DATA_CHUNKS:
        for col in order:
            data[col] = np.nan

        mapseq = mlt.datamap.MapSeq(drug_mode='smiles_1',
                                    protein_mode='protein_3',
                                    max_drug_len=200,
                                    max_protein_len=1000)

        drug_seqs = data.Drug.unique()
        protein_seqs = data.Target.unique()

        map_drug, map_protein = mapseq.create_maps(drug_seqs=drug_seqs, protein_seqs=protein_seqs)
        test_batch_size = mlt.training_utils.get_batch_size(data.shape[0], max_batch_size=MAX_BATCH_SIZE)
        test_gen = mlt.datagen.DataGen(data,
                                       map_drug,
                                       map_protein,
                                       shuffle=False,
                                       test_only=True)

        test_gen = test_gen.get_generator(test_batch_size)
        prediction = model.predict(test_gen, steps=data.shape[0] // test_batch_size, verbose=0)
        for k, col in enumerate(order[:-1]):
            data[col] = prediction[k]

        data.drop(['Drug', 'Target', 'qed'], axis=1).to_csv(path_result_save, mode='a', header=False)

    shutil.make_archive(path_result_save.replace('.csv', ''), 'zip', './', path_result_save)
    if os.path.exists(path_result_save):
        os.remove(path_result_save)


if __name__ == '__main__':
    inference()