"""
Functions to calculate the scores of matrices
"""

# Standard Library Imports
import pandas as pd
import numpy as np
# Local Library Imports

####################################################################################################
import math


def round_sig(x, sig=3):
    """local function to round a number to a given number of significant figures"""
    if x != 0:
        return round(x, sig-int(math.floor(math.log10(abs(x))))-1)
    else:
        return 0  # return 0 if the input is 0


def score_relative_matrix(sequence: str, matrix: pd.DataFrame) -> float:
    """Function that calculates the score of a sequence given, based on the relative matrix

    Args:
        matrix (pd.DataFrame): A pandas DataFrame with the relative matrix

    Returns:
        float: The score of the sequence given a matrix
    """

    # Initialize the score
    score = 1
    # For each position in the sequence
    for position in range(len(sequence)):
        # Get the nucleotide at this position
        nucleotide = sequence[position]
        # Get the probability of this nucleotide at this position
        prob = matrix.loc[position, nucleotide]
        # Multiply the score by the probability
        score *= prob
        # Round the score to 3 significant figures
    score = round_sig(score, 3)
    return score


def score_weight_matrix(sequence: str, matrix: pd.DataFrame) -> float:
    """Function that calculates the score of a sequence given, based on the weight matrix

    Args:
        matrix (pd.DataFrame): A pandas DataFrame with the weight matrix

    Returns:
        float: The score of the sequence given a matrix
    """

    # Initialize the score
    score = 0
    # For each position in the sequence
    for position in range(len(sequence)):
        # Get the nucleotide at this position
        nucleotide = sequence[position]
        # Get the weight of this nucleotide at this position
        weight = matrix.loc[position, nucleotide]
        # Multiply the score by the weight
        score += weight
        # Round the score to 3 significant figures
    score = round_sig(score, 3)
    return score


def max_score_relative_matrix(matrix: pd.DataFrame) -> float:
    """Function that calculates the maximum score of a matrix

    Args:
        matrix (pd.DataFrame): A pandas DataFrame with the relative matrix

    Returns:
        float: The maximum score of the matrix
    """
    # Get the maximum score of the matrix
    max_score = np.prod(matrix.max(axis=1))

    max_score = round_sig(max_score, 3)
    return max_score


def max_score_weight_matrix(matrix: pd.DataFrame) -> float:
    """Function that calculates the maximum score of a matrix

    Args:
        matrix (pd.DataFrame): A pandas DataFrame with the weight matrix

    Returns:
        float: The maximum score of the matrix
    """
    # Get the maximum score of the matrix
    max_score = matrix.max(axis=1).sum()
    # Round the score to 3 significant figures
    max_score = round_sig(max_score, 3)
    return max_score


def min_score_relative_matrix(matrix: pd.DataFrame) -> float:
    """Function that calculates the minimum score of a matrix

    Args:
        matrix (pd.DataFrame): A pandas DataFrame with the relative matrix

    Returns:
        float: The minimum score of the matrix
    """
    # Get the minimum score of the matrix
    min_score = np.prod(matrix.min(axis=1))
    # Round the score to 3 significant figures
    min_score = round_sig(min_score, 3)
    return min_score


def min_score_weight_matrix(matrix: pd.DataFrame) -> float:
    """Function that calculates the minimum score of a matrix

    Args:
        matrix (pd.DataFrame): A pandas DataFrame with the weight matrix

    Returns:
        float: The minimum score of the matrix
    """
    # Get the minimum score of the matrix
    min_score = matrix.min(axis=1).sum()
    # Round the score to 3 significant figures
    min_score = round_sig(min_score, 3)
    return min_score


def relative_score(min_score: float, max_score: float, score: float) -> float:
    """Function that calculates the relative score of a sequence given the minimum and maximum scores

    Args:
        min_score (float): The minimum score of the matrix
        max_score (float): The maximum score of the matrix
        score (float): The score of the sequence

    Returns:
        float: The relative score of the sequence
    """
    # Calculate the relative score
    relative_score = (score - min_score) / (max_score - min_score)
    # Round the score to 3 significant figures
    relative_score = round_sig(relative_score, 3)
    return relative_score
