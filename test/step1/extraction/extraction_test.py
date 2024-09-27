import os

from src.step1.extraction.mRNADataExtractor import mRNADataExtractor
from src.utils.utils import load_config

PROJECT_PATH = os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd())))

def main():
    config = load_config(PROJECT_PATH + '/config/config.yaml')
    data_extractor = mRNADataExtractor(config)
    data_extractor.run()

if __name__ == '__main__':
    main()
