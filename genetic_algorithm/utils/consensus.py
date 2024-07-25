import logging
from functools import lru_cache, wraps
from typing import Any, Callable, Iterable, Sequence
import numpy as np
from numba import njit, prange
from numpy.typing import NDArray
# from genetic_algorithm.config.log_config import CustomLogger

logger = logging.getLogger(__name__)
#### IMPORTANT ####
# All "dict types" of nucleotides follow this convention: (A, C, G, T) so any argument that requires
# a dict type of nucleotides should be passed as a tuple of 4 elements, for arrays and tuples
# Example: (1, 2, 3, 4) or (0.25, 0.25, 0.25, 0.25)

IUPACCODES = {
    (1, 0, 0, 0): 0,  # 'A'
    (0, 1, 0, 0): 1,  # 'C'
    (0, 0, 1, 0): 2,  # 'G'
    (0, 0, 0, 1): 3,  # 'T'
    (1, 1, 0, 0): 4,  # 'M'
    (1, 0, 1, 0): 5,  # 'R'
    (1, 0, 0, 1): 6,  # 'W'
    (0, 1, 1, 0): 7,  # 'S'
    (0, 1, 0, 1): 8,  # 'Y'
    (0, 0, 1, 1): 9,  # 'K'
    (1, 1, 1, 0): 10,  # 'V'
    (1, 1, 0, 1): 11,  # 'H'
    (1, 0, 1, 1): 12,  # 'D'
    (0, 1, 1, 1): 13,  # 'B'
    (1, 1, 1, 1): 14,  # 'N'
    (0, 0, 0, 0,): 14  # 'N
}

# Convert the IUPAC code dictionary to arrays for numba compatibility
IUPACKEYS = np.array([list(key)
                      for key in IUPACCODES], dtype=np.int32)
IUPACVALUES = np.array(
    [value for value in IUPACCODES.values()], dtype=np.int32)
IUPACCHARS = np.array(
    ['A', 'C', 'G', 'T', 'M', 'R', 'W', 'S', 'Y', 'K', 'V', 'H', 'D', 'B', 'N'], dtype='U1')


def print_helper(message):
    print(message)


def no_zeroes_decorator(func: Callable) -> Callable:
    @wraps(func)
    def wrapper(array: NDArray[np.float64], *args, **kwargs):
        no_zeroes = kwargs.pop('no_zeroes', True)
        if no_zeroes and isinstance(array, np.ndarray):
            array = np.where(array == 0, 0.8, array)
        return func(array, *args, **kwargs)
    return wrapper


def hash_decorator(func: Callable, sorted_surface: int = -1) -> Callable:
    """Decorator that ensures all inputs are hashable by converting iterables to tuples and sorting
    them based on sorted_surface.This only sorts args or kwargs that are sequences and have a depth 
    greater than or equal to sorted_surface. Applying this decoratorto a function that recieves
    heterogeneous arguments, or unsortable aruments will log an error, and not srt the argument.
    This decorator is useful for functions that require hashable arguments, such as caching 
    functions. Sorting is useful when the orderof the elements does not affect the result in the
    function.

    Args:
        func (Callable): Function to be decorated, that recieves hashable arguments.

        sorted_surface (int): The depth of the argument that start to be sorted. Default is -1,
                            none will be sorted.
    Returns: Callable: A wrapper function that converts its arguments to hashable and sorted tuples
                       before calling the original function.
    """

    @wraps(func)
    def wrapper(*args, **kwargs):
        hashable_args = convert_to_hashable_and_sort(args, 0, sorted_surface)
        hashable_kwargs = {k: convert_to_hashable_and_sort((v,), 0, sorted_surface)[
            0] for k, v in kwargs.items()}
        return func(*hashable_args, **hashable_kwargs)

    def convert_to_hashable_and_sort(args: Any, current_depth: int, sorted_surface: int):
        hashable_args: list[Any] = []
        for arg in args:
            if isinstance(arg, Iterable) and not isinstance(arg, (str, bytes, bytearray)):
                if current_depth >= sorted_surface >= 0 and isinstance(arg, Sequence):
                    try:
                        sorted_arg = sorted(
                            arg, key=lambda x: (str(type(x)), x))
                        hashable_args.append(tuple(convert_to_hashable_and_sort(
                            sorted_arg, current_depth + 1, sorted_surface)))
                    except TypeError:
                        logger.error("Could not sort the argument: %s", arg)
                        hashable_args.append(tuple(convert_to_hashable_and_sort(
                            arg, current_depth + 1, sorted_surface)))
                else:
                    # Convert without sorting
                    hashable_args.append(tuple(convert_to_hashable_and_sort(
                        arg, current_depth + 1, sorted_surface)))
            else:
                hashable_args.append(arg)
        return hashable_args

    return wrapper


def alignment_to_numpy_decorator(func: Callable) -> Callable:
    """
    Converts the alignment_array argument of a function to a numpy array, and the string to
    numerical values.This decorator takes a function that expects an alignment_array as a list of
    strings and modifies it to acceptthe alignment_array as a numpy array of strings. It is useful
    for functions that perform operations on sequencealignments where numpy array operations are
    more efficient.

    Args:
        func (Callable): Function to be decorated, that receives an alignment_array as a Sequence
                         of strings.

    Returns:
        Callable: A wrapper function that converts its alignment_array argument to a numpy array
                  before calling the original function.
    """
    @wraps(func)
    def wrapper(alignment_array: Sequence[str], *args, **kwargs):
        """
        Wrapper function that converts an alignment_array list to a numpy array.

        This function takes an alignment_array represented as a sequence of strings, converts it
        into a numpy array of numerical values, and then calls the original function with this array
        and any additional arguments.

        Args:
            alignment_array (Sequence[str]): The alignment_array to be converted, represented as a
                                            sequence of strings.
            *args: Variable length argument list for the decorated function.
            **kwargs: Arbitrary keyword arguments for the decorated function.

        Returns:
            The return value of the decorated function.
        """
        sequences_list = [list(sequence) for sequence in alignment_array]
        # Ensure each element is a single character
        sequences_array = np.array(sequences_list, dtype='U1')
        vectorized_ord = np.vectorize(ord)
        unicode_array: np.ndarray = vectorized_ord(sequences_array)
        return func(unicode_array, *args, **kwargs)
    return wrapper


@njit(parallel=True)
def nucleotide_counter_row(column: NDArray[np.int32]) -> NDArray[np.int32]:
    column_count = np.zeros(4, dtype=np.int32)
    # Count the number of each nucleotide in the column = (A, C, G, T)
    for nts in column:
        if nts not in [ord("A"), ord("C"), ord("G"), ord("T")]:
            raise ValueError(f"Invalid nucleotide found: {chr(nts)}")
        if nts == ord("A"):
            column_count[0] += 1
        elif nts == ord("C"):
            column_count[1] += 1
        elif nts == ord("G"):
            column_count[2] += 1
        elif nts == ord("T"):
            column_count[3] += 1
    return column_count


def absolute_matrix(alignment_array: NDArray[np.int32]) -> NDArray[np.int32]:
    """Function that generates an absolute matrix from a multiple sequence alignment_array

    Args:
        alignment_array (list): A list of sequences, previously aligned

    Returns:
        np.ndarray: A numpy array with the absolute matrix
    """
    # Initialize an empty array with columns for each nucleotide and dtype as int
    abs_matrix = np.zeros((len(alignment_array[0]), 4), dtype=np.int32)
    # For each position in the alignment_array
    for position in range(len(alignment_array[0])):
        # Get the nucleotide at this position in each sequence
        column = alignment_array[:, position]
        # Count the occurrences of each nucleotide
        abs_matrix[position] = nucleotide_counter_row(column)
    return abs_matrix


@njit
def _relativize_matrix(abs_matrix: NDArray[np.float64]) -> NDArray[np.float64]:
    """Core function that relativizes an absolute matrix using JIT compilation.

    Args:
        abs_matrix (np.ndarray): An absolute matrix

    Returns:
        np.ndarray: A relativized matrix
    """
    # Sum the values of each line
    row_sums = abs_matrix.sum(axis=1)
    row_sums_reshaped = row_sums.reshape(-1, 1)
    # Divide each value in the matrix by the sum of its line
    return np.round((abs_matrix / row_sums_reshaped), 2)


@no_zeroes_decorator
def relativize_matrix(abs_matrix: NDArray[np.int32]) -> NDArray[np.float64]:
    """Function that relativizes an absolute matrix

    Args:
        abs_matrix (np.ndarray): An absolute matrix

    Returns:
        np.ndarray: A relativized matrix
    """
    return _relativize_matrix(abs_matrix.astype(np.float64))


@njit
def calculate_weight(
        row: NDArray[np.float64], background: NDArray[np.float64]) -> NDArray[np.float64]:
    """Calculate the weight of a row in the matrix

    Args:
        row (np.ndarray): A row in the matrix
        background (np.ndarray): The background frequencies

    Returns:
        float: The weight of the row
    """
    for i in range(4):
        row[i] = row[i] * np.log2(row[i] / background[i])

    return row


def weight_matrix(rel_matrix: NDArray[np.int64], background: NDArray[np.float64]
                  = np.array((0.25, 0.25, 0.25, 0.25))) -> NDArray[np.float64]:

    weights = np.apply_along_axis(calculate_weight, 1, rel_matrix, background)
    return weights


def get_weight_matrix(alignment_array: NDArray[np.int32], background: NDArray[np.float64]
                      = np.array((0.25, 0.25, 0.25, 0.25))) -> NDArray[np.float64]:
    abs_matrix = absolute_matrix(alignment_array)
    rel_matrix = relativize_matrix(abs_matrix)
    return weight_matrix(rel_matrix, background)


def get_relative_matrix(alignment_array: NDArray[np.int32]) -> NDArray[np.float64]:
    abs_matrix = absolute_matrix(alignment_array)
    return relativize_matrix(abs_matrix)


@ njit(parallel=True)
def consensus(rel_matrix: NDArray[np.float64], threshold):
    sequence_length = rel_matrix.shape[0]
    consensus_sequence = np.empty(sequence_length, dtype=np.int32)
    for j in prange(sequence_length):  # pylint: disable=not-an-iterable
        top_freq = max(rel_matrix[j])
        selected = rel_matrix[j] / top_freq >= threshold
        key = (selected[0], selected[1], selected[2], selected[3])
        for k in range(len(IUPACKEYS)):  # pylint: disable=consider-using-enumerate
            if np.array_equal(IUPACKEYS[k], key):
                consensus_sequence[j] = IUPACVALUES[k]
                break
    return consensus_sequence


def convert_to_iupac_chars(consensus_sequence):
    return ''.join(IUPACCHARS[consensus_sequence])


@hash_decorator
@alignment_to_numpy_decorator
def alignment_to_consensus(alignment_array, threshold):
    rel_matrix = get_relative_matrix(alignment_array)
    consensus_sequence = consensus(rel_matrix, threshold)
    return convert_to_iupac_chars(consensus_sequence)
