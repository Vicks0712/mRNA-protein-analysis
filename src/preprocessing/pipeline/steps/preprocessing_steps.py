from src.logs.logger import logger


START_CODON = 'AUG'
STOP_CODONS = ['UAA', 'UAG', 'UGA']

# Step 1: Transcribir la secuencia
def transcribe_sequence(seq_record):
    seq_record.seq = seq_record.seq.transcribe()
    return seq_record

def find_start_codon(seq_record):
    seq_str = str(seq_record.seq)
    start_pos = seq_str.find(START_CODON)
    if start_pos == -1:
        logger.warning('PREPROCESSING-PIPELINE: Sequence %s without start codon', seq_record.id)
        return None
    seq_record.start_pos = start_pos
    return seq_record

# Step 3: Encontrar el codón de terminación
def find_stop_codon(seq_record):
    seq_str = str(seq_record.seq)
    stop_positions = [seq_str.find(codon, seq_record.start_pos) for codon in STOP_CODONS]
    stop_positions = [pos for pos in stop_positions if pos != -1]
    if not stop_positions:
        logger.warning('PREPROCESSING-PIPELINE: Sequence %s without end codon', seq_record.id)
        return None
    seq_record.stop_pos = min(stop_positions)
    return seq_record

# Step 4: Extraer la secuencia codificante
def extract_coding_sequence(seq_record):
    coding_seq = seq_record.seq[seq_record.start_pos:seq_record.stop_pos+3]
    seq_record.seq = coding_seq
    return seq_record