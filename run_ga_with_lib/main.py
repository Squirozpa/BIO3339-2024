
from dataclasses import dataclass
import random
import re
import Bio3339_tools as bt
from marboxes_benchmark import get_benchmark
from reader import read_fasta
import multiprocessing as mp

############################################ PARAMETERS ############################################

MARBOXES_GENES = read_fasta("_input_files/marboxes_27nts.fasta")
mutation_rate = 0.4
mutation_change = True
population_size = 100
immigration_rate = 0
gene_pool_type = bt.genetic_algorithm.strategies_classes.BinaryGenePool
survivor_strategy = bt.genetic_algorithm.strategies_classes.ElitismSurvivorSelection
cross_strategy = bt.genetic_algorithm.strategies_classes.TwoPointCrossover
parent_strategy = bt.genetic_algorithm.strategies_classes.RouletteWheelSelection
num_evolutions = 10
num_iterations = 1000
######################################### FITNESS FUNCTION #########################################

## Fitness Parameters ##

BENCHMARK = get_benchmark("_input_files/benchmark_marboxes.txt")
ALIGNMENT_DICT = {"A": ["A"], "C": ["C"], "T": ["T"], "G": ["G"], "N": ["A", "C", "G", "T"],
                  "R": ["A", "G"], "Y": ["C", "T"], "S": ["G", "C"],
                  "W": ["A", "T"], "K": ["G", "T"], "M": ["A", "C"], "B": ["C", "G", "T"],
                  "D": ["A", "G", "T"], "H": ["A", "C", "T"], "V": ["A", "C", "G"]}
SCORE_DICT = {"A": 2, "C": 2, "G": 2, "T": 2, "R": 1, "Y": 1, "S": 1, "W": 1, "K": 1, "M": 0,
              "B": 0, "D": 0, "H": 0, "V": 0, "N": -1}


@dataclass
class MarboxesFitness(bt.genetic_algorithm.FitnessStrategy):
    threshold: float = 0.6
    alignment_dict = ALIGNMENT_DICT
    score_dict = SCORE_DICT
    benchmark = BENCHMARK

    def _fitness(self, gene_vector: list) -> int:
        # Initialize score the base score if everything is N then score goes to 0
        score = 9000

        # Get the consensus sequence
        sequence_alignment = self.gene_pool.convert_genes(gene_vector)
        # print(sequence_alignment)
        consensus_sequence = bt.alignment_to_consensus(
            sequence_alignment, threshold=self.threshold)
        consensus_score = 0

        for nts in consensus_sequence:
            # calculates a score for the consensus sequence based on the score_dict
            consensus_score += self.score_dict[nts]
        # iterate over all the benchmarks
        for _, line in enumerate(self.benchmark):
            # iterate over the benchmark sequence (line[3])
            for i in range(len(line[3])-27):
                if i == line[2]:  # if the position is the same as the benchmark position
                    # if the alignment is not correct
                    if not self.alignment(line[3][i:i+27], consensus_sequence,):
                        # subrsact 27 from the score (one for every nucleotide)
                        score -= 27
                    else:
                        score += consensus_score  # if it aligns correctly add the "consensus score"
                else:
                    # if its not the benchmark position but aligns correctly
                    if self.alignment(line[3][i:i+27], consensus_sequence):
                        score -= 1  # subract 1 from the score
        return score

    def alignment(self, target, sequence):
        """Simple alignment function that returns True if the target sequence aligns with the sequence"""
        max_mismatch = 4
        for i in range(len(target)):
            if target[i] not in self.alignment_dict[sequence[i]]:
                max_mismatch -= 1
            if max_mismatch < 0:
                return False
        return True

######################################### GENETIC ALGORITHM ########################################


def initialize_sublasses(genes, seed):
    random_gen = random.Random(seed)
    gene_pool = gene_pool_type(genes=genes, random=random_gen)
    fitness_strategy = MarboxesFitness(gene_pool)
    crossover_strategy = cross_strategy(random=random_gen)
    survivor_selection_strategy = survivor_strategy(random=random_gen)
    parent_selection_strategy = parent_strategy(random=random_gen)
    return (
        gene_pool, fitness_strategy, crossover_strategy, survivor_selection_strategy,
        parent_selection_strategy, random_gen)


def initialize_population(genes, seed, shared_dict, population_state=None):
    (gene_pool, fitness_strategy, crossover_strategy, survivor_selection_strategy,
     parent_selection_strategy, random_gen) = initialize_sublasses(genes, seed)
    params = bt.genetic_algorithm.PopulationParams(mutation_rate=mutation_rate,
                                                   mutation_change=mutation_change,
                                                   population_size=population_size,
                                                   num_survivors=4, seed=seed,
                                                   random=random_gen)
    config = bt.genetic_algorithm.StrategyConfig(
        gene_pool, fitness_strategy, crossover_strategy, survivor_selection_strategy,
        parent_selection_strategy)
    if population_state is not None:
        return bt.Population(params=params, strategies=config, state=population_state, shared_dict=shared_dict)
    return bt.Population(params=params, strategies=config, state=None, shared_dict=shared_dict)


def initialize_population_state(individuals):
    return bt.genetic_algorithm.PopulationState(individuals=individuals)


def run_population(population: bt.Population, iterations: int):
    population.run(iterations)
    return population


def main(num_processes=3):
    manager = mp.Manager()
    shared_dict = manager.dict()

    with mp.Pool(num_processes) as pool:
        initial_populations = pool.starmap(initialize_population, [
            (MARBOXES_GENES, i, shared_dict) for i in range(num_processes)])
    print("All populations initialized")

    with mp.Pool(num_processes) as pool:
        results = pool.starmap(run_population, [
            (pop, num_evolutions) for pop in initial_populations])

    print("All populations ran")

    with mp.Pool(num_processes) as pool:
        results_states = pool.starmap(initialize_population_state, [
            ([results[i].state.individuals[:50] + results[i-1].state.individuals[:50]]) for i in range(len(results))])

    with mp.Pool(num_processes) as pool:
        results_new = pool.starmap(run_population, [
            (initialize_population(MARBOXES_GENES, 1, shared_dict, results_states[i]), num_evolutions) for i in range(len(results_states))])

    print("All populations ran")

    for result in results_new:
        print(result[0])

    for result in results:
        print(result[0])


if __name__ == "__main__":
    main()
