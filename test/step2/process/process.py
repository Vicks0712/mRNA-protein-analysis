import os

from src.step2.process.mRNAProcessorOrchester import mRNAProcessorOrchester
from src.utils.utils import load_config

PROJECT_PATH = os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd())))


def main(processing_key):
    """Main function to execute the mRNA translation process."""
    config = load_config(PROJECT_PATH + '/config/config.yaml')

    data_translator = mRNAProcessorOrchester(config, processing_key)
    data_translator.run()

if __name__ == '__main__':
    processing_key = 'protein_translation'
    main(processing_key)