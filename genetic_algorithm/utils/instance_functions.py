from logging import getLogger
import traceback
import json
from genetic_algorithm.abstract_classes.abstract_strategies import SurvivorSelectionStrategy
from genetic_algorithm.concrete_strategies.gene_pool_strategies import BinaryGenePool, IndexGenePool
from genetic_algorithm.concrete_strategies.parent_selection_strategies import RouletteWheelSelection
from genetic_algorithm.concrete_strategies.survivor_selection_strategies import ElitismSurvivorSelection
from genetic_algorithm.utils import reader
from genetic_algorithm.concrete_strategies.fitness_strategies import FitnessFunctionBinaryAmbiguity
from genetic_algorithm.concrete_strategies.fitness_strategies import FitnessFunctionBinaryBenchmark
from genetic_algorithm.concrete_strategies.crossover_strategies import TwoPointCrossover
from genetic_algorithm.population_class import Population, PopulationParams, PopulationState, StrategyConfig
from genetic_algorithm.utils import benchmark_marboxes
logger = getLogger(__name__)


def load_config(config_file: str) -> dict:
    """Load configuration from a JSON file."""
    config_file = 'genetic_algorithm/config/' + config_file
    try:
        with open(config_file, 'r', encoding='utf-8') as f:
            config = json.load(f)
        logger.debug("Configuration loaded from %s", config_file)
        return config
    except Exception as e:
        logger.error("Failed to load configuration: %s", e)
        raise


def init_population(
        name: str, shared_dict: dict[tuple, float],
        config_file: str, threshold: float = 0.1, **kwargs) -> Population:
    """Initialize a population with parameters from a config file."""
    config = load_config(config_file)

    # Initialize gene pool
    if not config['gene_pool_type'] == 'None':
        gene_pool = reader.read_fasta(config['gene_pool_file'])
    else:
        gene_pool = kwargs['gene_pool']
    gene_pool_strategy = globals()[config['gene_pool_type']](genes=gene_pool)

    # Create strategies based on config
    fitness_strategy = globals()[config['fitness_strategy']]
    survivor_selection_strategy = globals()[config['survivor_selection_strategy']]
    crossover_strategy = globals()[config['crossover_strategy']]
    parent_selection_strategy = globals()[config['parent_selection_strategy']]

    population_params = PopulationParams(
        mutation_rate=config['mutation_rate'],
        mutation_change=config['mutation_change'],
        population_size=config['population_size'],
        immigration_rate=config['immigration_rate'],
        seed=config['seed']
    )

    population_strategy = StrategyConfig(
        gene_pool=gene_pool_strategy,
        fitness=fitness_strategy,
        crossover=crossover_strategy,
        survivor_selection=survivor_selection_strategy,
        parent_selection=parent_selection_strategy,
        fitness_args={'threshold': threshold},
    )
    try:
        population = Population(
            params=population_params,
            strategies=population_strategy,
            state=None,
            name=name,
            shared_dict=shared_dict)
        logger.debug("Population initialized: %s", name)
        return population
    except Exception as e:
        logger.error("Error initializing population: %s", name)
        logger.error(e)
        traceback_str = traceback.format_exc()  # Get the full traceback as a string
        logger.error(traceback_str)  # Log the full traceback
        return None
