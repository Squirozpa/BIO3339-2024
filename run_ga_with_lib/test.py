
from operator import ge
import re
from main import initialize_sublasses, MarboxesFitness, initialize_population, run_population, MARBOXES_GENES
import Bio3339_tools as bt
import random
import multiprocessing as mp


def test_fitness():
    seeded_random = random.Random("dsadasda")
    gene_pool = bt.genetic_algorithm.strategies_classes.BinaryGenePool(
        MARBOXES_GENES, seeded_random)
    test_fitness = MarboxesFitness(gene_pool)
    print(test_fitness._fitness(
        [seeded_random.randint(0, 1) for _ in range(56)]))
    # Test initialize_subclasses


def test_initialize_subclasses(seed=1):
    genes = ["AC", "GC", "GC", "TC", "AC", "GC", "GC", "TC",
             "AC", "GC", "GC", "TC", "AC", "GC", "GC", "TC", "AC", "GC", "GC", "TC"]
    (gene_pool, fitness_strategy, crossover_strategy, survivor_selection_strategy,
     parent_selection_strategy, random_gen) = initialize_sublasses(
        genes, seed)
    assert isinstance(gene_pool, bt.genetic_algorithm.GenePool)
    assert isinstance(fitness_strategy, MarboxesFitness)
    assert isinstance(crossover_strategy,
                      bt.genetic_algorithm.CrossoverStrategy)
    assert isinstance(survivor_selection_strategy,
                      bt.genetic_algorithm.SurvivorSelectionStrategy)
    assert isinstance(parent_selection_strategy,
                      bt.genetic_algorithm.ParentSelectionStrategy)
    print("All tests passed!")
    random_list = []
    random_list.append(survivor_selection_strategy.random.randint(0, 100))
    random_list.append(gene_pool.random.randint(0, 100))
    random_list.append(parent_selection_strategy.random.randint(0, 100))
    return True, random_list


def test_initialize_population(seed=1, shared_dict=None):
    genes = ["AC", "GC", "GC", "TC", "AC", "GC", "GC", "TC",
             "AC", "GC", "GC", "TC", "AC", "GC", "GC", "TC", "AC", "GC", "GC", "TC"]
    population = initialize_population(
        genes, seed, shared_dict)
    assert isinstance(population, bt.Population)
    print("All tests passed!")
    return population


def test_run_genetic_algorithm(population):
    new_pop = run_population(population, 5)
    return new_pop


def main():
    manager = mp.Manager()
    shared_dict = manager.dict()

    with mp.Pool(3) as pool:
        results_1 = pool.map(test_initialize_subclasses, [1, 1, 1])
    for item in results_1:
        print(item[1])

    with mp.Pool(3) as pool:
        results_2 = pool.starmap(test_initialize_population, [
            (1, shared_dict) for _ in range(3)])
    for item in results_2:
        print(repr(item[1]))

    with mp.Pool(3) as pool:
        results_2 = pool.map(test_run_genetic_algorithm, results_2)
    for item in results_2:
        print(repr(item[0]))


if __name__ == "__main__":
    test_fitness()
