"""
Concrete Gene Pool Strategies for the Genetic Algorithm. These strategies are used to generate,
mutate and convert genes for the individuals in the population.
"""

# Standard Library Imports
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any
# Local Library Imports
from genetic_algorithm.abstract_classes.abstract_strategies import GenePoolStrategy
if TYPE_CHECKING:
    from genetic_algorithm.individual_class import Individual
####################################################################################################


@dataclass
class BinaryGenePool(GenePoolStrategy):
    genes: list
    vector_size: int = field(init=False)

    def __post_init__(self):
        self.vector_size = len(self.genes)

    def generate_genes(self) -> list:
        return [self.random.choice([0, 1]) for _ in range(self.vector_size)]

    def mutate_genes(self, individual: "Individual") -> "Individual":
        random_index = self.random.randint(0, self.vector_size - 1)
        individual.genes[random_index] = 1 - individual.genes[random_index]
        return individual

    def convert_genes(self, gene_vector: list) -> list:
        return [self.genes[i] for i in range(self.vector_size) if gene_vector[i] == 1]


@dataclass
class IndexGenePool(GenePoolStrategy):
    genes: list[Any]
    vector_size: int

    def generate_genes(self) -> list:
        genes = self.random.sample(
            range(len(self.genes)), self.vector_size)
        return genes

    def mutate_genes(self, individual: "Individual") -> "Individual":
        random_index = self.random.randint(0, self.vector_size - 1)
        possible_genes = [x for x in range(
            len(self.genes)) if x not in individual.genes]

        individual.genes[random_index] = self.random.choice(
            list(possible_genes))

        if len(individual.genes) != len(set(individual.genes)):
            print("Duplicate genes found in individual")

        return individual

    def convert_genes(self, gene_vector: list[int]) -> list:
        sorted_gene_vector = sorted(gene_vector)
        return [self.genes[i] for i in sorted_gene_vector]


@dataclass
class IndexChoiceGenePool(GenePoolStrategy):
    """Gene Pool Strategy, where each gene by index is a choice from a list of values. The gene is the index of the value in the gene pool."""
    genes: list[list[Any]]
    vector_size: int = field(init=False)

    def __post_init__(self):
        self.vector_size = len(self.genes)

    def generate_genes(self) -> list:
        gene = []
        for i in range(self.vector_size):
            gene.append(self.random.choice(self.genes[i]))
        return gene

    def mutate_genes(self, individual: "Individual") -> "Individual":
        random_index = self.random.randint(0, self.vector_size - 1)
        possible = set(self.genes[random_index]) - set(
            [individual.genes[random_index]])
        individual.genes[random_index] = self.random.choice(list(possible))
        return individual
