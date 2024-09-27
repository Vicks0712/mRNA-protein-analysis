from src.step2.process.mRNAProcessor import ARNmProcessor


class AtypicalTranslator(ARNmProcessor):
    def process_sequence(self, seq_record):
        # Ejemplo de traducción atípica, este método podría aplicar otro procesamiento.
        from Bio.Seq import Seq
        # A modo de ejemplo, puede realizar una "traducción" personalizada
        # basada en algún otro patrón o regla
        atypical_seq = Seq(seq_record.seq).reverse_complement().translate(to_stop=True)
        seq_record.seq = atypical_seq
        return seq_record