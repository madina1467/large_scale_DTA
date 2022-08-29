from sys import platform

USERNAME = 'madina_sultanova' #'maxim_popov'

if platform == "linux" or platform == "linux2": # linux
    HOME_PATH = '/home/' + USERNAME + '/large_scale_DTA/'
    DATA_PATH = '/raid/' + USERNAME + '/data/'
    DATA_PATH_BATCHED = DATA_PATH + 'batched/'
    DATA_PATH_PARQUET = DATA_PATH + 'batched_parquet/'
    DATA_PATH_PARQUET_PROCESSED = DATA_PATH_PARQUET + 'processed/'
    MODELS_PATH = '/raid/' + USERNAME + '/models/GraphDTA_Results/'
    RESULTS_PATH = DATA_PATH + 'results/'
    MLTLE_RESULTS_PATH = DATA_PATH + 'results_mltle/'
    FINAL_RESULTS_PATH = DATA_PATH + 'final_results/'
elif platform == "darwin": # OS X
    HOME_PATH = '/Users/madina/Documents/research/denovo_moleculardesign/Custom-GraphDTA-with-BindingDB/'
    DATA_PATH = HOME_PATH + 'large_scale_DTA/data/'
    DATA_PATH_BATCHED = DATA_PATH + 'batched/'
    DATA_PATH_PARQUET = DATA_PATH + 'batched_parquet/'
    DATA_PATH_PARQUET_PROCESSED = DATA_PATH_PARQUET + 'processed/'
    MODELS_PATH = HOME_PATH + 'GraphDTA_Results/'
    RESULTS_PATH = HOME_PATH + 'large_scale_DTA/results/'
    TEST_RESULTS_PATH_TEST = HOME_PATH + 'large_scale_DTA/test_results/'
    MLTLE_RESULTS_PATH = HOME_PATH + 'large_scale_DTA/real_res'
    FINAL_RESULTS_PATH = HOME_PATH + 'large_scale_DTA/final_results/'
elif platform == "win32": # Windows...
    pass

# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()