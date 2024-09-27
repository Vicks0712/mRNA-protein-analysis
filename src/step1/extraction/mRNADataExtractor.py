import os
from src.logs.logger import logger
from Bio import Entrez, SeqIO

class mRNADataExtractor:
    def __init__(self, config):
        self.project_path = config['project_path']
        self.email = config['email']
        self.query = config['data_extraction']['query']
        self.db = config['data_extraction']['db']
        self.rettype = config['data_extraction']['rettype']
        self.retmax = config['data_extraction']['retmax']
        self.output_dir = self._build_raw_data_path(config)
        Entrez.email = self.email

    def _build_raw_data_path(self, config):
        return self.project_path + config['data_extraction']['output_dir']

    def run(self):
        id_list = self._fetch_sequences()
        for seq_id in id_list:
            seq_record = self._download_sequence(seq_id)
            if seq_record:
                self._save_sequence(seq_record)

    def _fetch_sequences(self):
        logger.info('mRNA-DATA-EXTRACTOR: Initializing sequences extraction')
        handle = Entrez.esearch(db=self.db, term=self.query, retmax=self.retmax)
        record = Entrez.read(handle)
        handle.close()
        id_list = record["IdList"]
        logger.info('mRNA-DATA-EXTRACTOR: Founf %d sequences', len(id_list))
        return id_list

    def _download_sequence(self, seq_id):
        try:
            handle = Entrez.efetch(db=self.db, id=seq_id, rettype=self.rettype, retmode="text")
            seq_record = SeqIO.read(handle, "fasta")
            handle.close()
            return seq_record
        except Exception as e:
            logger.error('mRNA-DATA-EXTRACTOR: Error at downloading %s: %s', seq_id, str(e))
            return None

    def _save_sequence(self, seq_record):
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            logger.info('mRNA-DATA-EXTRACTOR: Directory created: %s', self.output_dir)
        filepath = os.path.join(self.output_dir, f"{seq_record.id}.fasta")
        SeqIO.write(seq_record, filepath, "fasta")
        logger.info('mRNA-DATA-EXTRACTOR: Sequence stored at %s', filepath)
