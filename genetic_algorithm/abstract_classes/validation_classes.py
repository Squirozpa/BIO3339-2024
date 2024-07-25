from random import Random
from typing import Optional
from genetic_algorithm.abstract_classes.abstract_strategies import (
    GenePoolStrategy, FitnessStrategy, CrossoverStrategy, SurvivorSelectionStrategy,
    ParentSelectionStrategy)


class ValidationStrategyConfig:

    @staticmethod
    def validate(gene_pool, fitness, crossover, survivor_selection, parent_selection):
        ValidationStrategyConfig._validate_gene_pool(gene_pool)
        ValidationStrategyConfig._validate_fitness(fitness)
        ValidationStrategyConfig._validate_crossover(crossover)
        ValidationStrategyConfig._validate_survivor_selection(survivor_selection)
        ValidationStrategyConfig._validate_parent_selection(parent_selection)

    @staticmethod
    def _validate_gene_pool(gene_pool: GenePoolStrategy):
        if not isinstance(gene_pool, GenePoolStrategy):
            raise ValueError(
                f"gene_pool must be an instance of GenePoolStrategy. Received {gene_pool}")

    @staticmethod
    def _validate_fitness(fitness):
        if not isinstance(fitness, FitnessStrategy):
            raise ValueError(
                f"fitness must be an instance of FitnessStrategy. Received {fitness}")

    @staticmethod
    def _validate_crossover(crossover: CrossoverStrategy):
        if not isinstance(crossover, CrossoverStrategy):
            raise ValueError(
                f"crossover must be an instance of CrossoverStrategy. Received {crossover}")

    @staticmethod
    def _validate_survivor_selection(survivor_selection: SurvivorSelectionStrategy):
        if not isinstance(survivor_selection, SurvivorSelectionStrategy):
            raise ValueError(
                f"survivor_selection must be an instance of SurvivorSelectionStrategy. "
                f"Received {survivor_selection}")

    @staticmethod
    def _validate_parent_selection(parent_selection: ParentSelectionStrategy):
        if not isinstance(parent_selection, ParentSelectionStrategy):
            raise ValueError(
                f"parent_selection must be an instance of ParentSelectionStrategy. "
                f"Received {parent_selection}")


class ValidationPopulationParams:

    @staticmethod
    def validate(
            mutation_rate, mutation_change, population_size,
            immigration_rate, num_survivors, seed, random):
        ValidationPopulationParams._validate_mutation_rate(mutation_rate)
        ValidationPopulationParams._validate_mutation_change(mutation_change)
        ValidationPopulationParams._validate_population_size(population_size)
        ValidationPopulationParams._validate_immmigration_rate(immigration_rate)
        ValidationPopulationParams._validate_num_survivors(num_survivors)
        ValidationPopulationParams._validate_seed(seed)
        ValidationPopulationParams._validate_random(random)

    @staticmethod
    def _validate_mutation_rate(mutation_rate: float):
        if not isinstance(mutation_rate, float):
            raise ValueError(
                f"mutation_rate must be a float. Received {mutation_rate}")
        if mutation_rate < 0 or mutation_rate > 1:
            raise ValueError(
                f"mutation_rate must be between 0 and 1. Received {mutation_rate}")

    @staticmethod
    def _validate_mutation_change(mutation_change: bool | None):
        if not isinstance(mutation_change, bool):
            if mutation_change is not None:
                raise ValueError(
                    f"mutation_change must be a boolean or None. Received {mutation_change}")

    @staticmethod
    def _validate_population_size(population_size: int):
        if not isinstance(population_size, int):
            raise ValueError(
                f"population_size must be an integer. Received {population_size}")

    @staticmethod
    def _validate_immmigration_rate(immigration_rate: Optional[float]):
        if immigration_rate is not None:
            if not isinstance(immigration_rate, float):
                raise ValueError(f"immigration_rate must be a float. Received{immigration_rate}")

    @staticmethod
    def _validate_num_survivors(num_survivors: int):
        if not isinstance(num_survivors, int):
            raise ValueError(f"num_survivors must be an integer. Received {num_survivors}")
        if num_survivors < 0:
            raise ValueError(
                f"num_survivors must be greater or equal than 0. Received {num_survivors}")

    @staticmethod
    def _validate_seed(seed: Optional[int | str | float]):
        if seed is not None:
            if not isinstance(seed, (int, str, float)):
                raise ValueError(f"seed must be an int, str, or float. Received {seed}")

    @staticmethod
    def _validate_random(random: Optional[Random]):
        if random is not None:
            if not isinstance(random, Random):
                raise ValueError(f"random must be an instance of Random. Received {random}")
