"""
Concrete Crossover Strategies for the Genetic Algorithm. These strategies are used to crossover the
genes of two parents to create new children.
"""

# Standard Library Imports
from ctypes import cast
from dataclasses import dataclass

from random import Random
from typing import Iterable, Sequence
# Local Library Imports
from genetic_algorithm.abstract_classes.abstract_strategies import CrossoverStrategy
from genetic_algorithm.individual_class import Individual

####################################################################################################


@dataclass
class SinglePointCrossover(CrossoverStrategy):
    # Default Single Point Crossover Method returns two children from two parents
    def crossover(self, parents: Iterable[Individual]) -> Iterable[Individual]:
        parent1, parent2 = tuple(parents)
        children = []
        crossover_point = self.random.randint(0, len(parent1.genes))
        child1_genes = parent1.genes[:crossover_point] + \
            parent2.genes[crossover_point:]
        child2_genes = parent2.genes[:crossover_point] + \
            parent1.genes[crossover_point:]
        child1 = Individual(genes=child1_genes, iteration=self.population.state.iteration,
                            population_id=self.population.id)
        child2 = Individual(genes=child2_genes, iteration=self.population.state.iteration,
                            population_id=self.population.id)
        children.append(child1)
        children.append(child2)
        return children


@dataclass
class TwoPointCrossover(CrossoverStrategy):
    # Default Two Point Crossover Method returns two children from two parents
    def crossover(self, parents: Iterable[Individual]):
        parent1, parent2 = tuple(parents)
        children = []
        crossover_points = sorted([self.random.randint(
            0, len(parent1.genes)), self.random.randint(0, len(parent1.genes))])

        child1_genes = parent1.genes[:crossover_points[0]] + parent2.genes[crossover_points[0]:crossover_points[1]] + parent1.genes[crossover_points[1]:]
        child2_genes = parent2.genes[:crossover_points[0]] + parent1.genes[crossover_points[0]:crossover_points[1]] + parent2.genes[crossover_points[1]:]
        child1 = Individual(genes=child1_genes, iteration=self.population.state.iteration,
                            population_id=self.population.id)
        child2 = Individual(genes=child2_genes, iteration=self.population.state.iteration,
                            population_id=self.population.id)
        children.append(child1)
        children.append(child2)
        return children


@dataclass
class TwoPointCrossOverForUniqueIndex(CrossoverStrategy):
    def crossover(self, parents: Sequence[Individual]) -> Iterable[Individual]:
        if len(parents) != 2:
            raise ValueError(
                "TwoPointCrossOverForUniqueIndex requires exactly 2 parents")
        parent1_genes, parent2_genes = tuple(
            parent.genes for parent in parents)
        # Extract the common genes between the two parents
        common_genes = set(parent1_genes).intersection(set(parent2_genes))
        unique_genes_parent1 = list(set(parent1_genes) - common_genes)
        unique_genes_parent2 = list(set(parent2_genes) - common_genes)

        # If there are no unique genes, return the parents
        if not unique_genes_parent1 or not unique_genes_parent2:
            return parents

        crossover_points = sorted([self.random.randint(0, len(
            unique_genes_parent1) - 1), self.random.randint(0, len(unique_genes_parent2) - 1)])

        child1_genes = unique_genes_parent1[:crossover_points[0]] + unique_genes_parent2[crossover_points[0]:crossover_points[1]] + unique_genes_parent1[crossover_points[1]:]
        child2_genes = unique_genes_parent2[:crossover_points[0]] + unique_genes_parent1[crossover_points[0]:crossover_points[1]] + unique_genes_parent2[crossover_points[1]:]
        child1_genes.extend(common_genes)
        child2_genes.extend(common_genes)
        return Individual(genes=child1_genes, iteration=self.population.state.iteration, population_id=self.population.id), Individual(genes=child2_genes, iteration=self.population.state.iteration, population_id=self.population.id)


@dataclass
class UniformCrossover(CrossoverStrategy):
    # Default Uniform Crossover Method returns two children from two parents
    def crossover(self, parents: Sequence[Individual]):
        parent1, parent2 = tuple(parents)
        children = []
        child1_genes = []
        child2_genes = []
        for gene1, gene2 in zip(parent1.genes, parent2.genes):
            if self.random.choice([True, False]):
                child1_genes.append(gene1)
                child2_genes.append(gene2)
            else:
                child1_genes.append(gene2)
                child2_genes.append(gene1)
        child1 = Individual(genes=child1_genes, iteration=self.population.state.iteration,
                            population_id=self.population.id)
        child2 = Individual(genes=child2_genes, iteration=self.population.state.iteration,
                            population_id=self.population.id)
        children.append(child1)
        children.append(child2)
        return children
