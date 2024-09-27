from abc import abstractmethod, ABC


class ARNmProcessor(ABC):
    @abstractmethod
    def process_sequence(self, seq_record):
        pass