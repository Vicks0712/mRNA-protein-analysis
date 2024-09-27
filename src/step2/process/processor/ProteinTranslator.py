from src.step2.process.mRNAProcessor import ARNmProcessor
from Bio.Seq import Seq


class ProteinTranslator(ARNmProcessor):
    def process_sequence(self, seq_record):
        protein_seq = Seq(seq_record.seq).translate(to_stop=True)
        seq_record.seq = protein_seq
        return seq_record