"""
Class for population of individuals. This class is designed to hold a list of individuals, which its
genes are binary vectors. Also contains all the logic to generate a new population from the current 
one, such as selection, crossover, and mutation. The class is designed to be used in a genetic
algorithm, but can be used in any optimization algorithm that requires a binary representation of 
the solution. It contains all the logic to run the genetic algorithm, refer to the README for more 
information on how to use it.
"""

# Standard Library Imports

from logging import getLogger
from random import Random
from typing import Optional, cast
# Local Library Imports
from genetic_algorithm.utils import consensus
from genetic_algorithm.config.log_config import CustomLogger
from genetic_algorithm.abstract_classes.abstract_strategies import (
    CrossoverStrategy, ParentSelectionStrategy, SurvivorSelectionStrategy,
    FitnessStrategy)
from genetic_algorithm.individual_class import Individual
from genetic_algorithm.abstract_classes.abstract_population import (
    AbstractPopulation, StrategyConfig, PopulationState, PopulationParams)
####################################################################################################
logger = getLogger(__name__)


class Population(AbstractPopulation):
    def __init__(self,
                 params: PopulationParams,
                 strategies: StrategyConfig,
                 state: Optional[PopulationState],
                 name: Optional[str] = None,
                 shared_dict: Optional[dict] = None
                 ):
        logger.debug("Initializing concrete population: %s", name)
        super().__init__(strategies, params, state, name, shared_dict)
        self._generate_initial_population()
        logger.debug("Concrete population initialized: %s", name)

    def __str__(self):
        return (f"Population ID: {self.id}, Size: {self.params.population_size}, Best Individual: "
                f"{self.state.best_individual}, Iterations: {self.state.iteration}, "
                f" Seed: {self.params.seed}, Total Fitness: {self.state.total_fitness}")

    def _generate_initial_population(self) -> None:
        """Generates the initial population of individuals."""
        for _ in range(self.params.population_size):
            individual = self._generate_random_individual()
            self.append(individual)
        self._compute_fitness()
        self.state.refresh_state()

    def _compute_fitness(self):
        """Computes the fitness of all the individuals in the population."""
        self.strategies.fitness = cast(FitnessStrategy, self.strategies.fitness)
        for individual in self:
            individual.fitness = self.strategies.fitness.fitness(individual.genes)

    def _evolve_population(self):
        """Evolves the population to the next generation."""
        self.strategies.survivor_selection = cast(
            SurvivorSelectionStrategy, self.strategies.survivor_selection)
        self.params.random = cast(Random, self.params.random)
        self.strategies.crossover = cast(CrossoverStrategy, self.strategies.crossover)
        self.strategies.parent_selection = cast(
            ParentSelectionStrategy, self.strategies.parent_selection)
        self.strategies.parent_selection = cast(
            ParentSelectionStrategy, self.strategies.parent_selection)

        survivors = list(self.strategies.survivor_selection.select_survivors(
            self.params.num_survivors))
        new_population = survivors

        while len(new_population) < self.params.population_size:
            parents = self.strategies.parent_selection.select_parents(self)
            children = self.strategies.crossover(parents)

            for child in children:
                if self.params.random.random() < self.params.mutation_rate:
                    child = self.strategies.gene_pool.mutate_genes(child)
                new_population.append(child)

        self.state.individuals = new_population[:self.params.population_size]
        self._compute_fitness()
        self.sort()
        self.state.refresh_state()
        self.state.add_iteration()

    def _run(self, iterations: int):
        self._evolve_population()
        if self.params.mutation_change:
            self.params.mutation_rate = self._mutation_rate_change(iterations)
        print(f"Population {self.name} has finished running an iteration.")

    def _mutation_rate_change(self, iterations: int):
        """Changes the mutation rate based on the iteration."""
        self.logger = cast(CustomLogger, self.logger)

        if self.params.mutation_rate < 0.01:
            self.logger.debug("Mutation rate is already at its minimum.")
            self.params.mutation_change = False
            return 0.01
        new_mutation_rate = self.params.mutation_rate * (1 - self.state.iteration / iterations)
        self.logger.debug("Mutation rate has been changed to: %s", str(new_mutation_rate))
        return new_mutation_rate

    def _generate_random_individual(self):
        """Generates a random individual."""
        gene_vector = self.strategies.gene_pool.generate_genes()
        return Individual(genes=gene_vector, population_id=self.id, iteration=self.state.iteration)
