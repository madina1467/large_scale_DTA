import csv

import pandas as pd
from rdkit import Chem
from tqdm.auto import tqdm
import os
from sys import platform
# if platform == "darwin": # OS X
#     from large_scale_DTA.const import *
# else: # linux
from const import *

from Bio import SeqIO


def create_dataset():

    if not os.path.exists(DATA_PATH_BATCHED):
        os.makedirs(DATA_PATH_BATCHED)
        print("The new directory is created: ", DATA_PATH_BATCHED)
    if not os.path.exists(DATA_PATH_PARQUET):
        os.makedirs(DATA_PATH_PARQUET)
        print("The new directory is created: ",DATA_PATH_PARQUET)

    ligands = csv.reader(open(DATA_PATH + "drugbank_db_smiles.csv"))
    fasta_sequences = SeqIO.parse(open(
        DATA_PATH + "uniprot-compressed_true_download_true_format_fasta_query__2A_20AND_2-2022.07.01-07.13.49.42.fasta"),
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

    if platform == "darwin": # OS X
        n_prots = 20  # 20386   20
        n_drugs = 12  # 10920   12
    else: # linux
        n_prots = len(prots) #20386   20
        n_drugs = len(drugs) #10920   12

    row_pairs = n_prots * n_drugs

    iterations = 48
    batch_size = row_pairs // iterations

    printProgressBar(0, iterations, prefix='Progress:', suffix='Complete', length=iterations)

    count, prot_ind, drug_ind = 0, 0, 0

    for iter in range(iterations): #tqdm(range(iterations), desc = f'Progress:', position=0, leave=True, ascii=True, ncols= iterations*2):
        breaker = False
        dataset = str(iter) + '_drugbank_uniprot_' + str(batch_size)

        with open(DATA_PATH_BATCHED + dataset + '.csv', 'w') as f:
            f.write('Drug,Target,Y\n')  # Drug,Target,Y

            pbar = tqdm(total=batch_size, position=0,  leave=True, ascii=True, )
            while prot_ind < n_prots:
                while drug_ind < n_drugs:
                    count = count + 1
                    ls = []
                    ls += [drugs[drug_ind]]
                    ls += [prots[prot_ind]]
                    ls += [1]
                    f.write(','.join(map(str, ls)) + '\n')

                    pbar.update(1)
                    # print(count, dataset, prot_ind, drug_ind)

                    if drug_ind == n_drugs - 1:
                        drug_ind = 0
                        prot_ind = prot_ind + 1
                    else:
                        drug_ind = drug_ind + 1

                    if count == batch_size:
                        count = 0
                        breaker = True
                        break
                if breaker:
                    break
            pbar.close()

        df = pd.read_csv(DATA_PATH_BATCHED + dataset + '.csv')
        df.to_parquet(DATA_PATH_PARQUET + dataset + '.parquet')

        printProgressBar(iter + 1, iterations, prefix='Progress:', suffix='Complete', length=iterations)

if __name__ == '__main__':
    create_dataset()
