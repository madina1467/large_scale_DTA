import glob
import os
import pandas as pd

from fastparquet import ParquetFile


from large_scale_DTA.const import RESULTS_PATH, FINAL_RESULTS_PATH, TEST_RESULTS_PATH_TEST


def transform_result_files():
    if not os.path.exists(FINAL_RESULTS_PATH):
        os.makedirs(FINAL_RESULTS_PATH)
        print("The new directory is created: ", FINAL_RESULTS_PATH)

    parquet_files = glob.glob(RESULTS_PATH + '*.parquet')
    for file in parquet_files:
        filename = os.path.splitext(os.path.basename(file))[0]
        print('Transform file: ', filename)

        # df = pd.read_parquet(file)
        # df.to_csv(FINAL_RESULTS_PATH + filename + '.csv')

        pf = ParquetFile(file)
        dataFrame = pf.to_pandas()
        dataFrame.to_csv(FINAL_RESULTS_PATH + filename + '.csv', index=False)

def shrink_real_results():
    if not os.path.exists(TEST_RESULTS_PATH_TEST):
        os.makedirs(TEST_RESULTS_PATH_TEST)
        print("The new directory is created: ", TEST_RESULTS_PATH_TEST)

    parquet_files = glob.glob(TEST_RESULTS_PATH_TEST + '*.parquet')
    for file in parquet_files:
        filename = os.path.splitext(os.path.basename(file))[0]
        print('Transform file: ', filename)

        # df = pd.read_parquet(file)
        # df.to_csv(FINAL_RESULTS_PATH + filename + '.csv')

        pf = ParquetFile(file)
        dataFrame = pf.to_pandas()

        dataFrame = dataFrame[0:5]
        # dataFrame.to_csv(FINAL_RESULTS_PATH + filename + '.csv', index=False)

        dataFrame.to_parquet(TEST_RESULTS_PATH_TEST + str(filename) + '_' + 'res.parquet', index=False)


def create_new_data():
    if not os.path.exists(TEST_RESULTS_PATH_TEST):
        os.makedirs(TEST_RESULTS_PATH_TEST)
        print("The new directory is created: ", TEST_RESULTS_PATH_TEST)

    parquet_files = glob.glob(TEST_RESULTS_PATH_TEST + '*.parquet')
    for file in parquet_files:
        filename = os.path.splitext(os.path.basename(file))[0]
        print('Transform file: ', filename)

        pf = ParquetFile(file)
        dataFrame = pf.to_pandas()

        col_names = list(dataFrame.columns)

        for col in col_names:
            if col.lower() != 'drug' and col.lower() != 'target' and col.lower() != 'gcnnet_bdtdc_kd':
                dataFrame.drop(col, axis=1, inplace=True)

        dataFrame['test'] = 'test_val'

        # dataFrame.to_csv(FINAL_RESULTS_PATH + filename + '.csv', index=False)

        dataFrame.to_parquet(TEST_RESULTS_PATH_TEST + str(filename) + '_' + 'fake_add.parquet', index=False)

if __name__ == '__main__':
    create_new_data()



