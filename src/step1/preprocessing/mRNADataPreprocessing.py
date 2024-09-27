
import os
from Bio import SeqIO
from src.step1.preprocessing.pipeline.preprocessing_pipeline import PreprocessingPipeline
from src.logs.logger import logger

class mRNADataPreprocessing:

    _project_path = os.path.dirname(os.path.dirname(os.getcwd()))

    def __init__(self, config):
        self.project_path = config['project_path']
        self.input_dir = self._build_data_path(config['data_cleaning']['input_dir'])
        self.output_dir = self._build_data_path(config['data_cleaning']['output_dir'])
        self.pipeline = PreprocessingPipeline()

    def _build_data_path(self, data_path):
        return self.project_path + data_path

    def run(self):
        self.validate_output_directory()
        self.process_sequences()
        logger.info('mRNA-DATA-PREPROCESSING: Sequence cleaning process completed')

    def validate_output_directory(self):
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            logger.info('mRNA-DATA-PREPROCESSING: Created output directory: %s', self.output_dir)
        else:
            logger.info('mRNA-DATA-PREPROCESSING: Output directory already exists: %s', self.output_dir)

    def process_sequences(self):
        for filename in os.listdir(self.input_dir):
            filepath = os.path.join(self.input_dir, filename)
            self.process_file(filepath)

    def process_file(self, filepath):
        try:
            seq_record = self.read_sequence(filepath)
            if seq_record:
                cleaned_seq_record = self.pipeline.run(seq_record)  # Run the pipeline on the sequence
                self.log_processing_result(filepath, seq_record, cleaned_seq_record)
                if cleaned_seq_record:
                    self.save_cleaned_sequence(cleaned_seq_record)
        except Exception as e:
            logger.error(f'mRNA-DATA-PREPROCESSING: Error processing file {filepath}: {str(e)}')

    def read_sequence(self, filepath):
        try:
            seq_record = SeqIO.read(filepath, "fasta")
            logger.info(f'mRNA-DATA-PREPROCESSING: Reading sequence {seq_record.id} from {filepath}')
            return seq_record
        except Exception as e:
            logger.error(f'mRNA-DATA-PREPROCESSING: Failed to read sequence from {filepath}: {str(e)}')
            return None

    def save_cleaned_sequence(self, cleaned_seq_record):
        output_path = os.path.join(self.output_dir, f"{cleaned_seq_record.id}_clean.fasta")
        SeqIO.write(cleaned_seq_record, output_path, "fasta")
        logger.info(f'mRNA-DATA-PREPROCESSING: Saved cleaned sequence: {output_path}')

    def log_processing_result(self, filepath, seq_record, cleaned_seq_record):
        if cleaned_seq_record:
            logger.info(f'mRNA-DATA-PREPROCESSING: Successfully processed sequence {seq_record.id} from {filepath}')
        else:
            logger.warning(f'mRNA-DATA-PREPROCESSING: Sequence {seq_record.id} from {filepath} could not be cleaned')
