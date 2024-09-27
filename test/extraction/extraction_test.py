import os
import yaml

from src.extraction.mRNADataExtractor import mRNADataExtractor

PROJECT_PATH = os.path.dirname(os.getcwd())

def load_config(config_file):
    with open(config_file, 'r') as file:
        config = yaml.safe_load(file)
    return config

def main():
    config = load_config(PROJECT_PATH + '/config/config.yaml')
    data_extractor = mRNADataExtractor(config)
    data_extractor.run()

if __name__ == '__main__':
    main()
