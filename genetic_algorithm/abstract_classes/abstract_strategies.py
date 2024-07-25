from __future__ import annotations
from abc import ABC, abstractmethod

from multiprocessing.managers import DictProxy
from typing import Any, Iterable, TYPE_CHECKING, Sequence
from random import Random

from dataclasses import dataclass, field

if TYPE_CHECKING:
    from genetic_algorithm.abstract_classes.abstract_population import AbstractPopulation
    from genetic_algorithm.individual_class import Individual


@dataclass
class GenePoolStrategy(ABC):
    """This class contains the info of a gene pool. It is important that the random_gene function is
    consistent with the FitnessStrategy as its the core functionality of this algorithm"""
    genes: list = field(repr=False)
    vector_size: int = field()
    random: Random = field(default=Random(), repr=False, compare=False)

    def __len__(self):
        return self.vector_size

    def __getitem__(self, index):
        return self.genes[index]

    @abstractmethod
    def generate_genes(self) -> list:
        pass

    @abstractmethod
    def mutate_genes(self, individual: "Individual") -> "Individual":
        pass

    @abstractmethod
    def convert_genes(self, gene_vector: list) -> Iterable:
        pass


@dataclass
class FitnessStrategy(ABC):
    gene_pool: GenePoolStrategy
    shared_dict: dict[tuple, float] | None = None

    def fitness(self, gene_vector) -> float:
        fitness_value = self._fitness(gene_vector)
        if self.shared_dict is not None:
            self.shared_dict[tuple(gene_vector)] = fitness_value
        return fitness_value

    @abstractmethod
    def _fitness(self, gene_vector) -> float:
        pass


@dataclass
class CrossoverStrategy(ABC):
    population: AbstractPopulation | None = None
    random: Random | None = Random()
    args: tuple = field(default_factory=tuple)
    kwargs: dict = field(default_factory=dict)

    def __call__(self, *args: Any, **kwds: Any) -> Any:
        return self.crossover(*args, **kwds)

    @abstractmethod
    def crossover(self, parents: Sequence['Individual']) -> Iterable[Individual]:
        pass


@dataclass
class SurvivorSelectionStrategy(ABC):
    population: AbstractPopulation | None = None
    random: Random | None = Random()
    args: tuple = field(default_factory=tuple)
    kwargs: dict = field(default_factory=dict)

    @abstractmethod
    def select_survivors(self, num_survivors) -> Iterable[Individual]:
        pass


@dataclass
class ParentSelectionStrategy(ABC):
    population: AbstractPopulation | None = None
    random: Random = Random()
    args: tuple = field(default_factory=tuple)
    kwargs: dict = field(default_factory=dict)

    @abstractmethod
    def select_parents(self, *args, **kwargs) -> Sequence[Individual]:
        pass
