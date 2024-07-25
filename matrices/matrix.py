"""
Functions to generate the matrices from a given multiple sequence alignment
"""

# Standard Library Imports
import pandas as pd
import numpy as np
# Local Library Imports

####################################################################################################


def reverse_complement(sequence: str) -> str:
    """Internal function te generate the reverse complement of a DNA sequence."""
    swap_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join([swap_dict[char] for char in sequence[::-1]])


def genomic_count(sequence: str) -> dict:
    genomic_count = {"A": 0, "C": 0, "G": 0, "T": 0}
    reversed_sequence = reverse_complement(sequence)
    for nucleotide in sequence:
        genomic_count[nucleotide] += 1
    for nucleotide in reversed_sequence:
        genomic_count[nucleotide] += 1
    genomic_frequency = {key: value / (2 *
                         len(sequence)) for key, value in genomic_count.items()}
    return genomic_frequency


def absolute_matrix(alignment: tuple) -> pd.DataFrame:
    """Function that generates a absolute matrix from a multiple sequence alignment

    Args:
        alignment (list): A list of sequences, previosly aligned

    Returns:
        pd.DataFrame: A pandas DataFrame with the absolute matrix
    """

    # Create an empty DataFrame with columns for each nucleotide
    abs_matrix = pd.DataFrame(columns=['A', 'C', 'G', 'T'])
    # For each position in the alignment
    for position in range(len(alignment[0])):
        # Get the nucleotide at this position in each sequence
        nucleotides = [seq[position] for seq in alignment]
        # Count the occurrences of each nucleotide
        counts = pd.Series(nucleotides).value_counts()
        # Append the counts to the absolute matrix
        abs_matrix = abs_matrix._append(counts, ignore_index=True)
    # Fill any missing values with 0
    abs_matrix = abs_matrix.fillna(0)

    return abs_matrix


def relative_matrix(alignment: tuple, no_zeroes: bool = False) -> pd.DataFrame:
    """Function that generates a relative matrix from a multiple sequence alignment

    Args:
        alignment (list): A list of sequences, previosly aligned
        no_zeroes (bool, optional): If the relative matrix should have zeroes. Defaults to False. Used for calculating the weight matrix.
    Returns:
        pd.DataFrame: A pandas DataFrame with the relative matrix
    """

    # Get the absolute matrix
    abs_matrix = absolute_matrix(alignment)
    # Get the total number of sequences
    total_seqs = len(alignment)
    # If the relative matrix should have zeroes
    if no_zeroes:
        abs_matrix[abs_matrix == 0] = 0.8
    # Divide each element in the absolute matrix by the total number of sequences
    rel_matrix = abs_matrix / total_seqs
    # Convert each element in the relative matrix to float
    rel_matrix = rel_matrix.astype(float)
    # Round each element in the relative matrix to 2 decimals
    rel_matrix = rel_matrix.round(2)
    return rel_matrix


def entropy_matrix(alignment: list) -> pd.DataFrame:
    """Function that generates a entropy matrix from a multiple sequence alignment

    Args:
        alignment (list): A list of sequences, previosly aligned

    Returns:
        pd.DataFrame: A pandas DataFrame with the entropy matrix
    """
    # Get the relative matrix
    rel_matrix = relative_matrix(alignment, no_zeroes=True)
    # Create a copy of the relative matrix
    entropy_matrix = rel_matrix.copy()
    # For each column in the relative matrix
    for column in rel_matrix.columns:
        # Calculate the entropy of the column
        entropy = -1 * (rel_matrix[column] * rel_matrix[column].apply(
            lambda x: 0 if x == 0 else np.log2(x))).sum()
        # Calculate the entropy of the column
        entropy = 2 - entropy
        # Multiply each element in the column by the entropy
        entropy_matrix[column] = rel_matrix[column] * entropy
    # Round each element in the entropy matrix to 2 decimals
    entropy_matrix = entropy_matrix[column].round(2)
    return entropy_matrix


def weight_matrix(alignment: tuple, genome_frequencey: dict = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}) -> pd.DataFrame:
    """Function that generates a weight matrix from a multiple sequence alignment

    Args:
        alignment (list): A list of sequences, previosly aligned

    Returns:
        pd.DataFrame: A pandas DataFrame with the weight matrix
    """
    # Get the relative matrix
    rel_matrix = relative_matrix(alignment, no_zeroes=True)
    # Create a copy of the relative matrix
    weight_matrix = rel_matrix.copy()
    # For each column in the relative matrix
    for column in rel_matrix.columns:
        # Calculate the weight of the column
        weight_matrix[column] = rel_matrix[column] / genome_frequencey[column]
        weight_matrix[column] = weight_matrix[column].apply(
            lambda x: np.log2(x))
    weight_matrix = weight_matrix.round(2)
    return weight_matrix
