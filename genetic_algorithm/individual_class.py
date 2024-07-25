"""
Contains the class for each individual in the population. The individual's genes are represented by a binary vector and has a fitness value associated with it. The class also contains multiple methods to interact with the individual, such as comparison, iteration, and representation. The class is designed
to be used in a genetic algorithm, but can be used in any optimization algorithm that requires a binary representation of the solution. It is only a class for holding the information of the individual, and does not contain any optimization algorithm logic.
"""

# Standard Library Imports
import copy
from dataclasses import dataclass, field
import logging
from typing import Optional, Any

# Local Library Imports

####################################################################################################


@dataclass
class Individual:
    """Class for Individual in a Genetic Algorithm."""
    iteration: int
    genes: list[Any]
    population_id: int
    _id: int = field(default=0, init=False)
    _fitness: Optional[float | int] = field(
        default=None, compare=False, init=False)

    @property
    def fitness(self) -> float | int:
        """Fitness getter method."""
        if self._fitness is None:
            raise ValueError("Fitness value is None")
        return self._fitness

    @fitness.setter
    def fitness(self, value: float):
        """Fitness setter method."""
        if not isinstance(value, (int, float)):
            raise ValueError("Fitness value must be a number")
        self._fitness = value

    @property
    def id(self) -> int:
        """Id getter method."""
        return self._id

    def __post_init__(self):
        self.logger = logging.getLogger(__name__)
        self._id = id(self)
        self._population_id = self.population_id

    # region Magic Methods
    def __str__(self):
        return f'Individual: {self.id}, Fitness: {self.fitness}, Iteration {self.iteration}'

    def __repr__(self):
        return f'Individual: {self.id}, Fitness: {self.fitness}, Iteration {self.iteration},' + \
            f'Population: {self.population_id}, Genes: {self.genes}'

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Individual):
            return NotImplemented
        return self.genes == other.genes

    @staticmethod
    def _validate_other(other: "Individual"):
        if not isinstance(other, Individual):
            raise ValueError("Can't compare individual with non-individual")
        if other.fitness is None:
            raise ValueError(
                "Can't compare individuals with None fitness values")
        return other.fitness

    def __lt__(self, other: "Individual") -> bool:
        other_fitness = self._validate_other(other)

        return self.fitness < other_fitness

    def __le__(self, other: "Individual") -> bool:
        other_fitness = self._validate_other(other)

        if self.fitness == other_fitness:
            return self.iteration >= other.iteration
        else:
            return self.fitness < other_fitness

    def __gt__(self, other: "Individual"):
        other_fitness = self._validate_other(other)
        return self.fitness > other_fitness

    def __ge__(self, other: "Individual"):
        other_fitness = self._validate_other(other)
        if self.fitness == other.fitness:
            return self.iteration <= other.iteration
        return self.fitness > other_fitness

    def __ne__(self, other: object):
        if not isinstance(other, Individual):
            raise ValueError("Can't compare individual with non-individual")
        return self.genes != other.genes

    def __len__(self):
        return len(self.genes)

    def __getitem__(self, key):
        return self.genes[key]

    def __setitem__(self, key, value):
        self.genes[key] = value

    def __iter__(self):
        return iter(self.genes)

    def __contains__(self, item):
        return item in self.genes

    def __add__(self, other: "Individual"):
        if isinstance(other, Individual):
            if self.fitness is not None and other.fitness is not None:
                return self.fitness + other.fitness
            else:
                raise ValueError(
                    "Can't add individuals with None fitness values")

    def __mul__(self, other: "Individual"):
        if isinstance(other, Individual):
            if self.fitness is not None and other.fitness is not None:
                return self.fitness * other.fitness
            else:
                raise ValueError(
                    "Can't multiply individuals with None fitness values")
        else:
            raise ValueError(
                "Can't multiply individual with non-individual")

    def __rmul__(self, other: "Individual"):
        if isinstance(other, Individual):
            if self.fitness is not None and other.fitness is not None:
                return self.fitness * other.fitness
            else:
                raise ValueError(
                    "Can't multiply individuals with None fitness values")

    def __iadd__(self, other: "Individual"):
        if isinstance(other, Individual):
            if self.fitness is not None and other.fitness is not None:
                self.fitness += other.fitness
                return self
            else:
                raise ValueError(
                    "Can't add individuals with None fitness values")

    def deepcopy(self, memo=None):
        """Deepcopy method."""
        new_indidivual = copy.deepcopy(self, memo)  # type: ignore
        new_indidivual._id = id(  # pylint: disable=protected-access
            new_indidivual)
        return new_indidivual

    def __bool__(self):
        return bool(self.genes)

    def __int__(self):
        return int(''.join(str(bit) for bit in self.genes), 2)

    # endregion
