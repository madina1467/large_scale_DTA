import glob
import os
import pandas as pd

from fastparquet import ParquetFile


from large_scale_DTA.const import RESULTS_PATH, FINAL_RESULTS_PATH


def trnsform_result_files():
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


if __name__ == '__main__':
    trnsform_result_files()