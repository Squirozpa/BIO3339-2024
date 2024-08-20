"""
Main module to run the genetic algorithm.
"""
# Standard Library Imports
import argparse
import logging.handlers
import multiprocessing
from typing import Optional

from numpy import sort
# Local Library Imports
from genetic_algorithm.concrete_strategies.gene_pool_strategies import BinaryGenePool, IndexGenePool
from genetic_algorithm.concrete_strategies.fitness_strategies import (
    FitnessFunctionBinaryAmbiguity, FitnessFunctionIndexAmbiguity)
from genetic_algorithm.concrete_strategies.crossover_strategies import (
    TwoPointCrossOverForUniqueIndex, TwoPointCrossover)
from genetic_algorithm.concrete_strategies.parent_selection_strategies import RouletteWheelSelection
from genetic_algorithm.concrete_strategies.survivor_selection_strategies import (
    ElitismSurvivorSelection)
from genetic_algorithm.config.log_config import setup_logging
from genetic_algorithm.population_class import (
    Population, StrategyConfig, PopulationParams, PopulationState)
from genetic_algorithm.utils import reader, consensus
from genetic_algorithm.utils.instance_functions import init_population
from genetic_algorithm.utils.worker_functions import run_population
from genetic_algorithm.utils import consensus
####################################################################################################

####################################################################################################


def main(config_file: str, num_population=2, max_workers: Optional[int] = 2, timeout: int = 300, runs: int = 10):
    """Main function to run the genetic algorithm."""
    logger = logging.getLogger(__name__)
    logger.info("Starting the genetic algorithm.")
    logger.debug("Starting the multiprocessing manager")

    manager = multiprocessing.Manager()
    shared_dict = manager.dict()

    logger.debug("Initializing the populations.")
    with multiprocessing.Pool(max_workers) as pool:
        init_populations = pool.starmap(init_population, [
            (f"Population_{i}", shared_dict, config_file, i / 10 + 0.1) for i in range(num_population)])
    logger.debug("Completed population initialization.")

    logger.debug("Running the populations.")
    with multiprocessing.Pool(max_workers) as pool:
        results = pool.starmap(run_population, [
            (population, runs) for population in init_populations])

    results = sorted(
        results, key=lambda x: x.state.best_individual.fitness, reverse=True)  # type: ignore
    for result in results:
        logger.info(result)
        logger.info(result.state.best_individual)
        logger.info(result.state.best_individual.genes)

    with open("_output_files/results.txt", "w") as file:
        with open("_input_files/marboxes_27nts.fasta", "r") as fasta_file:
            fasta = fasta_file.readlines()
            genes_counter = [[fasta[i], 0] for i in range(0, len(fasta), 2)]
            fasta_file.close()
        for result in results:
            file.write(f"{result}\n")
            file.write(f"{result.state.best_individual}\n")
            file.write(f"{consensus.alignment_to_consensus(result.strategies.gene_pool.convert_genes(result.state.best_individual.genes), 0.6)}\n")
            for i, gene in enumerate(result.state.best_individual.genes):
                if gene == 1:
                    genes_counter[i][1] += 1
        file.write("Gene count \n")
        for gene in genes_counter:

            file.write(f"{gene[0]}: {gene[1]}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the genetic algorithm.")
    parser.add_argument("config_file", type=str, help="Configuration file to use.")
    parser.add_argument(
        "--num_population", type=int, default=2, help="Number of populations to run.")
    parser.add_argument("--max_workers", type=int, default=2, help="Number of workers to use.")
    parser.add_argument("--timeout", type=int, default=300, help="Timeout for the workers.")
    parser.add_argument("--runs", type=int, default=10, help="Number of iterations to run.")
    args = parser.parse_args()

    main(config_file=args.config_file,
         num_population=args.num_population,
         max_workers=args.max_workers,
         timeout=args.timeout,
         runs=args.runs)
