class ClassMapper:
    _map = {
        "protein_translation": ProteinTranslator,
        "atypical_translation": AtypicalTranslator,
        "regulatory_processing": RegulatoryProcessor,
    }

    @staticmethod
    def get_class(key: str) -> type:
        return ClassMapper._map.get(key, None)