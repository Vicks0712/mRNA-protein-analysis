import os

from src.step1.preprocessing.mRNADataPreprocessing import mRNADataPreprocessing
from src.utils.utils import load_config

PROJECT_PATH = os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd())))


def main():
    config = load_config(PROJECT_PATH + '/config/config.yaml')

    data_preprocessor = mRNADataPreprocessing(config)
    data_preprocessor.run()

if __name__ == '__main__':
    main()
