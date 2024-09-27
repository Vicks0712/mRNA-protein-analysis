from src.step2.process.processor.AtypicalTranslator import AtypicalTranslator
from src.step2.process.processor.ProteinTranslator import ProteinTranslator
from src.step2.process.processor.RegulatoryProcessor import RegulatoryProcessor


class ProcessorMapper:
    _map = {
        "protein_translation": ProteinTranslator,
        "atypical_translation": AtypicalTranslator,
        "regulatory_processing": RegulatoryProcessor,
    }

    @staticmethod
    def get_class(key: str) -> type:
        return ProcessorMapper._map.get(key, None)