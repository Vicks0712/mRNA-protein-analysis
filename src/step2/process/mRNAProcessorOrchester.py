import os
from Bio import SeqIO
from src.step2.process.factory.ProcessorFactory import ProcessorFactory
from src.logs.logger import logger


class mRNAProcessorOrchester:
    def __init__(self, config, processing_key):
        # Verifica el directorio del proyecto
        self.project_path = config['project_path']

        # Configura los directorios de entrada y salida
        self.input_dir = self._build_data_path(config['sequence_translation']['input_dir'])
        self.output_dir = self._build_data_path(config['sequence_translation']['output_dir'])

        # Inicializa el tipo de procesamiento
        self.processing_key = processing_key
        self.processor = ProcessorFactory.initialize_class(self.processing_key)

        # Agregar logging para verificar los directorios
        logger.info(f'Project path: {self.project_path}')
        logger.info(f'Input directory: {self.input_dir}')
        logger.info(f'Output directory: {self.output_dir}')

    def _build_data_path(self, data_path):
        # Combina el path del proyecto con el data_path
        return self.project_path + data_path

    def run(self):
        """Inicia el proceso de traducción de secuencias."""
        self.validate_output_directory()
        self.process_sequences()
        logger.info('mRNA-DATA-TRANSLATION: Sequence translation process completed')

    def validate_output_directory(self):
        """Crea el directorio de salida si no existe."""
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            logger.info(f'mRNA-DATA-TRANSLATION: Created output directory: {self.output_dir}')
        else:
            logger.info(f'mRNA-DATA-TRANSLATION: Output directory already exists: {self.output_dir}')

    def process_sequences(self):
        """Procesa todas las secuencias en el directorio de entrada."""
        logger.info(f'mRNA-DATA-TRANSLATION: Starting to process files from {self.input_dir}')
        for filename in os.listdir(self.input_dir):
            filepath = os.path.join(self.input_dir, filename)
            logger.info(f'mRNA-DATA-TRANSLATION: Processing file {filename}')
            self.process_file(filepath)

    def process_file(self, filepath):
        """Procesa un archivo de secuencia individual."""
        try:
            seq_record = self.read_sequence(filepath)
            if seq_record:
                translated_seq_record = self.processor.process_sequence(seq_record)
                self.log_processing_result(filepath, seq_record, translated_seq_record)
                if translated_seq_record:
                    self.save_translated_sequence(translated_seq_record)
        except Exception as e:
            logger.error(f'mRNA-DATA-TRANSLATION: Error processing file {filepath}: {str(e)}')

    def read_sequence(self, filepath):
        """Lee un archivo de secuencia."""
        try:
            seq_record = SeqIO.read(filepath, "fasta")
            logger.info(f'mRNA-DATA-TRANSLATION: Reading sequence {seq_record.id} from {filepath}')
            return seq_record
        except Exception as e:
            print(f"PATH {filepath}")
            logger.error(f'mRNA-DATA-TRANSLATION: Failed to read sequence from {filepath}: {str(e)}')
            return None

    def save_translated_sequence(self, translated_seq_record):
        """Guarda la secuencia traducida en el directorio de salida."""
        output_path = os.path.join(self.output_dir, f"{translated_seq_record.id}_translated.fasta")
        SeqIO.write(translated_seq_record, output_path, "fasta")
        logger.info(f'mRNA-DATA-TRANSLATION: Saved translated sequence: {output_path}')

    def log_processing_result(self, filepath, seq_record, translated_seq_record):
        """Registra el resultado del proceso de traducción."""
        if translated_seq_record:
            logger.info(f'mRNA-DATA-TRANSLATION: Successfully translated sequence {seq_record.id} from {filepath}')
        else:
            logger.warning(f'mRNA-DATA-TRANSLATION: Sequence {seq_record.id} from {filepath} could not be translated')
