import multiprocessing
import logging
from genetic_algorithm.config.log_config import setup_logging


def simple_worker(num):
    logger = logging.getLogger(__name__)
    logger.debug(f"Worker: {num}")
    return num


def simple_test():
    logger = logging.getLogger(__name__)
    logger.info("Starting simple test")

    with multiprocessing.Pool(2) as pool:
        results = pool.map(simple_worker, range(5))
    logger.info(f"Results: {results}")


if __name__ == "__main__":
    simple_test()
