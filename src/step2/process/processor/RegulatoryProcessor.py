from src.step2.process.mRNAProcessor import ARNmProcessor


class RegulatoryProcessor(ARNmProcessor):
    def process_sequence(self, seq_record):
        # Simulación de un proceso regulatorio en lugar de traducción a proteínas
        # Solo como ejemplo, podría simplemente retornar una parte del ARN
        regulatory_seq = seq_record.seq[:30]  # Se podría realizar un análisis de regulación
        seq_record.seq = regulatory_seq
        return seq_record