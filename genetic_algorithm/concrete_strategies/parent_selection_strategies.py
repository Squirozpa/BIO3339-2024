"""
Selection Strategies are used to select parents from the population to be used in the crossover process.
"""

# Standard Library Imports
from dataclasses import dataclass
from typing import cast
from typing import TYPE_CHECKING, Sequence
# Local Library Imports
from genetic_algorithm.abstract_classes.abstract_strategies import ParentSelectionStrategy
if TYPE_CHECKING:
    from genetic_algorithm.individual_class import Individual
    from genetic_algorithm.population_class import Population
####################################################################################################


@dataclass
class RouletteWheelSelection(ParentSelectionStrategy):
    def select_parents(self, *args, **kwargs) -> Sequence["Individual"]:
        self.population = cast("Population", self.population)
        parents = []
        for _ in range(2):  # Repeat the selection process twice
            total_fitness = self.population.state.total_fitness
            r = self.random.uniform(0, total_fitness)
            s = 0.0
            parent_selected = False
            for individual in self.population:
                s += individual.fitness
                if s >= r:
                    parents.append(individual)
                    parent_selected = True
                    break
# Fallback in case no individual is selected (should not happen if population is correctly set up)
            if not parent_selected:
                parents.append(self.population.state.individuals[-1])
        return parents


@dataclass
class TournamentSelection(ParentSelectionStrategy):
    tournament_size: int = 6

    def select_parents(self, *args, **kwargs) -> Sequence["Individual"]:
        self.population = cast("Population", self.population)
        parents = []
        for _ in range(2):
            tournament = self.random.sample(
                self.population.state.individuals, self.tournament_size)
            parents.append(max(individual for individual in tournament))
        return parents
