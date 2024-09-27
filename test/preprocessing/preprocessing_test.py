import yaml
import os

from src.step1.preprocessing import mRNADataPreprocessing

PROJECT_PATH = os.path.dirname(os.path.dirname(os.getcwd()))

def load_config(config_file):
    with open(config_file, 'r') as file:
        config = yaml.safe_load(file)
    return config

def main():
    config = load_config(PROJECT_PATH + '/config/config.yaml')

    data_preprocessor = mRNADataPreprocessing(config)
    data_preprocessor.run()

if __name__ == '__main__':
    main()
