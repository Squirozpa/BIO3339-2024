"""
Concrete Fitness Strategies, here the fitness function should be defined based on the problem to
optimize, still this implementation works for the current problem, but to use it in other problems,
the fitness function should be changed.
"""
# Standard Library Imports
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Sequence
# Local Library Imports
from genetic_algorithm.abstract_classes.abstract_strategies import FitnessStrategy
from genetic_algorithm.utils import consensus
from genetic_algorithm.utils.benchmark_marboxes import get_benchmark
if TYPE_CHECKING:
    from genetic_algorithm.individual_class import Individual
####################################################################################################


@dataclass
class FitnessFunctionBinaryAmbiguity(FitnessStrategy):
    threshold: float = field(default=0.6)

    def _fitness(self, gene_vector) -> float:
        """
        Fitness function for the ambiguity problem using binary genes
        """
        alignment = self.gene_pool.convert_genes(gene_vector)
        consensus_sequence = consensus.alignment_to_consensus(
            alignment, threshold=self.threshold)
        unique = 3
        double = 2
        triple = 1
        quad = 0
        score = 0
        score_dic = {"A": unique, "C": unique, "T": unique, "G": unique, "R": double,
                     "Y": double, "S": double, "W": double, "K": double, "M": triple, "B": triple,
                     "D": triple, "H": triple, "V": triple, "N": quad}
        for nts in consensus_sequence:
            if nts not in score_dic:
                raise ValueError(f"Key '{nts}' not found in score_dic.")
            score += score_dic[nts]
        return score


@dataclass
class FitnessFunctionIndexAmbiguity(FitnessStrategy):
    threshold: float = 0.6
    prioritize_upper: bool = True

    def fitness(self, gene_vector) -> float:
        """
        Fitness function for the ambiguity problem using index genes
        """
        alignment = self.gene_pool.convert_genes(gene_vector)
        consensus_sequence = consensus.alignment_to_consensus(
            alignment, threshold=self.threshold)
        unique = 3
        double = 2
        triple = 1
        quad = 0
        score = 0
        score_dic = {"A": unique, "C": unique, "T": unique, "G": unique, "R": double,
                     "Y": double, "S": double, "W": double, "K": double, "M": triple,
                     "B": triple, "D": triple, "H": triple, "V": triple, "N": quad}

        for nts in consensus_sequence:
            if nts not in score_dic:
                raise ValueError(f"Key '{nts}' not found in score_dic.")
            score += score_dic[nts]
        return score


@dataclass
class FitnessFunctionBinaryBenchmark(FitnessStrategy):
    threshold: float = 0.6
    benchmark = get_benchmark("_input_files/benchmark_marboxes.txt")
    alignment_dict = {"A": ["A"], "C": ["C"], "T": ["T"], "G": ["G"], "N": ["A", "C", "G", "T"],
                      "R": ["A", "G"], "Y": ["C", "T"], "S": ["G", "C"],
                      "W": ["A", "T"], "K": ["G", "T"], "M": ["A", "C"], "B": ["C", "G", "T"],
                      "D": ["A", "G", "T"], "H": ["A", "C", "T"], "V": ["A", "C", "G"]}
    score_dict = {"A": 2, "C": 2, "G": 2, "T": 2, "R": 1, "Y": 1, "S": 1, "W": 1, "K": 1, "M": 0,
                  "B": 0, "D": 0, "H": 0, "V": 0, "N": -1}

    def _fitness(self, gene_vector) -> int:
        """
        Fitness function for the Benchmark problem using binary genes
        """
        score = 9000
        if self.benchmark is None:
            raise ValueError("Benchmark not provided.")

        alignment = self.gene_pool.convert_genes(gene_vector)
        consensus_sequence = consensus.alignment_to_consensus(
            alignment, threshold=self.threshold)
        consensus_score = 0
        for nt in consensus_sequence:
            consensus_score += self.score_dict[nt]
        for _, line in enumerate(self.benchmark):
            for j in range(len(line[3]) - 27):
                if j == line[2]:
                    if not self.alignment(line[3][j:j + 27], consensus_sequence):
                        score -= 27
                    else:
                        score += consensus_score
                else:
                    if self.alignment(consensus_sequence, line[3][j:j + 27]):
                        score -= 1
        return score

    def alignment(self, target, sequence):
        max_mismatch = 2
        for i in range(len(target)):
            if target[i] not in self.alignment_dict[sequence[i]]:
                max_mismatch -= 1
            if max_mismatch < 0:
                return False
        return True


@dataclass
class FitnessKnapsack(FitnessStrategy):
    knapsack_capacity: int = 15

    def fitness(self, gene_vector: Sequence[Sequence[int]]) -> float:
        """
        Fitness function for the Knapsack problem
        """
        weight = 0
        value = 0
        for gene in gene_vector:
            weight += gene[0]
            value += gene[1]
        if weight > self.knapsack_capacity:
            return 0
        return value
