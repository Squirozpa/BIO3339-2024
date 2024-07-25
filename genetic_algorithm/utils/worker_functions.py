from multiprocessing.managers import BaseManager
from genetic_algorithm.population_class import Population


def run_population(population: Population, iterations: int):
    """Runs the genetic algorithm for a population."""
    population.run(iterations)
    return population
