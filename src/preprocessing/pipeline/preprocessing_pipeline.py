from src.preprocessing.pipeline.steps.preprocessing_steps import transcribe_sequence, find_start_codon, find_stop_codon, extract_coding_sequence


class PreprocessingPipeline:
    def __init__(self):
        self.steps = [
            transcribe_sequence,
            find_start_codon,
            find_stop_codon,
            extract_coding_sequence
        ]

    def run(self, seq_record):
        for step in self.steps:
            seq_record = step(seq_record)
            if seq_record is None:
                break
        return seq_record