import logging

# Define custom logging levels
RESULTS = 27
BENCHMARK = 25
VERBOSE = 15
TOTAL_BENCHMARK = 7

logging.addLevelName(RESULTS, "RESULTS")
logging.addLevelName(BENCHMARK, "BENCHMARK")
logging.addLevelName(VERBOSE, "VERBOSE")
logging.addLevelName(TOTAL_BENCHMARK, "TOTAL_BENCHMARK")

# Custom logger class


class CustomLogger(logging.Logger):
    def __init__(self, name, level=logging.NOTSET):
        super().__init__(name, level)

    def results(self, message, *args, **kwargs):
        self.log(RESULTS, message, *args, **kwargs)

    def benchmark(self, message, *args, **kwargs):
        self.log(BENCHMARK, message, *args, **kwargs)

    def verbose(self, message, *args, **kwargs):
        self.log(VERBOSE, message, *args, **kwargs)

    def total_benchmark(self, message, *args, **kwargs):
        self.log(TOTAL_BENCHMARK, message, *args, **kwargs)


# Register the custom logger class
logging.setLoggerClass(CustomLogger)

# Configure logging


def setup_logging(value: int = logging.DEBUG):
    # Create custom loggers
    logger = logging.getLogger()
    logger.setLevel(value)

    # Create console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)
    console_formatter = logging.Formatter('%(message)s')
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

    # Create file handlers
    results_handler = logging.FileHandler('result_population.log')
    results_handler.setLevel(RESULTS)
    results_formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s')
    results_handler.setFormatter(results_formatter)

    benchmark_handler = logging.FileHandler('benchmark.log')
    benchmark_handler.setLevel(BENCHMARK)
    benchmark_formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s')
    benchmark_handler.setFormatter(benchmark_formatter)

    total_benchmark_handler = logging.FileHandler('total_benchmark.log')
    total_benchmark_handler.setLevel(TOTAL_BENCHMARK)
    total_benchmark_formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s')
    total_benchmark_handler.setFormatter(total_benchmark_formatter)

    # Add handlers to logger
    logger.addHandler(results_handler)
    logger.addHandler(benchmark_handler)
    logger.addHandler(total_benchmark_handler)

    numba_logger = logging.getLogger('numba')
    numba_logger.setLevel(logging.WARNING)


# Call setup_logging to configure the loggers
setup_logging(TOTAL_BENCHMARK)
