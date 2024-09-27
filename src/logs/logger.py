import logging
import os

class Logger:
    _instance = None

    def __new__(cls, log_file: str = 'mRNA_analysis.log', log_level: int = logging.DEBUG):
        if cls._instance is None:
            cls._instance = super(Logger, cls).__new__(cls)
            cls._instance._initialize(log_file, log_level)
        return cls._instance

    def _initialize(self, log_file: str, log_level: int):
        """
        Initialize the logger by creating file and console handlers, setting format, and attaching handlers.
        """
        self.logger = self._create_logger(log_level)
        formatter = self._create_formatter()
        file_handler = self._create_file_handler(log_file, formatter)
        console_handler = self._create_console_handler(formatter)

        self._attach_handlers(file_handler, console_handler)

    def _create_logger(self, log_level: int) -> logging.Logger:
        """
        Create a logger instance with the specified log level.

        Args:
            log_level (int): The log level for the logger.

        Returns:
            logging.Logger: The created logger.
        """
        logger = logging.getLogger("hermes-video")
        logger.setLevel(log_level)
        return logger

    def _create_formatter(self) -> logging.Formatter:
        """
        Create a log formatter.

        Returns:
            logging.Formatter: The log formatter.
        """
        return logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    def _create_file_handler(self, log_file: str, formatter: logging.Formatter) -> logging.FileHandler:
        """
        Create a file handler to log messages to a file.

        Args:
            log_file (str): The log file name.
            formatter (logging.Formatter): The log formatter.

        Returns:
            logging.FileHandler: The file handler for logging to a file.
        """
        log_directory = self._create_log_directory()
        file_handler = logging.FileHandler(os.path.join(log_directory, log_file))
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        return file_handler

    def _create_console_handler(self, formatter: logging.Formatter) -> logging.StreamHandler:
        """
        Create a console handler to log messages to the console.

        Args:
            formatter (logging.Formatter): The log formatter.

        Returns:
            logging.StreamHandler: The console handler for logging to the console.
        """
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        console_handler.setFormatter(formatter)
        return console_handler

    def _create_log_directory(self) -> str:
        """
        Create a directory for log files if it doesn't exist.

        Returns:
            str: The path to the log directory.
        """
        log_directory = os.path.dirname(os.getcwd()) + '/logs'  # Corregir ruta del directorio de logs
        if not os.path.exists(log_directory):
            os.makedirs(log_directory)
        return log_directory

    def _attach_handlers(self, file_handler: logging.FileHandler, console_handler: logging.StreamHandler):
        """
        Attach file and console handlers to the logger.

        Args:
            file_handler (logging.FileHandler): The file handler.
            console_handler (logging.StreamHandler): The console handler.
        """
        if not self.logger.hasHandlers():
            self.logger.addHandler(file_handler)
            self.logger.addHandler(console_handler)

    def get_logger(self):
        return self.logger


logger = Logger().get_logger()
