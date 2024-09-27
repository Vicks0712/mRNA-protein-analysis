from src.step2.process.mRNAProcessor import ARNmProcessor
from src.step2.process.factory.ProcessorMapper import ProcessorMapper


class ProcessorFactory:
    @staticmethod
    def initialize_class(key: str) -> ARNmProcessor:
        processor_class = ProcessorMapper.get_class(key)
        if processor_class:
            return processor_class()
        else:
            raise ValueError(f"Class not found for the given key: {key}")