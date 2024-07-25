"""
Functions to generate consensus sequences from weight or relative matrices
"""

# Standard Library Imports
import pandas as pd
import numpy as np
# Local Library Imports

####################################################################################################


def iupac_ambiguity(nucleotides: list) -> str:
    """Internal funtion that returns the IUPAC ambiguity code for a list of nucleotides"""
    iupac_dic = {"AG": "R", "CT": "Y", "CG": "S",
                 "AT": "W", "GT": "K", "AC": "M", "CGT": "B", "AGT": "D", "ACT": "H", "ACG": "V"}

    nucleotides_str = "".join(sorted(nucleotides))
    return iupac_dic.get(nucleotides_str)


def consensus(matrix: pd.DataFrame, matrix_type: str = "weight", threshold: float = 0.5, prioritize_upper: bool = False) -> str:
    """Function that generates a consensus sequence from a weight or relative matrix

    Args:
        matrix (pd.DataFrame): A pandas DataFrame with the weight or relative matrix
        threshold (int, optional): The threshold to consider a nucleotide in the consensus sequence
    Returns:
        str: The consensus sequence
    """
    # Thresholds are calculated based on the "ideal" relative frequency of nucleotides, ej 1, for single, 0.5 for double, 0.33 for triple
    # Then that frequency is divided by 0.25 then log2 is applied to the result, this is the maximum weight or relative for "ideal" nucleotides
    # The "ideal" weight or relative is then multiplied by the threshold to get the threshold for the weight or relative matrix
    # All above each thresholds are considered for that condition.
    if matrix_type == "weight":
        lower_threshold = np.log2(1.32) * threshold
        middle_threshold = threshold
        upper_threshold = 2 * threshold
    elif matrix_type == "relative":
        lower_threshold = 0.25 + (threshold*0.08)
        middle_threshold = 0.25 + (threshold*0.25)
        upper_threshold = 0.25 + (threshold*0.75)
    # Initialize the consensus sequence
    consensus = ''
    # For each position in the matrix
    for position in range(len(matrix)):
        # Finds nucleotides that exceed the threshold
        nucleotides = pd.to_numeric(matrix.iloc[position], errors='coerce')
        nucleotides_above_lower = nucleotides[nucleotides >= lower_threshold]
        nucleotides_above_middle = nucleotides[nucleotides >= middle_threshold]
        nucleotides_above_upper = nucleotides[nucleotides >= upper_threshold]
        # Starts appending with the upper threshold
        if prioritize_upper:
            if len(nucleotides_above_upper) == 1:
                consensus += nucleotides_above_upper.index[0]
            elif len(nucleotides_above_middle) == 2:
                consensus += iupac_ambiguity(
                    nucleotides_above_middle.index.tolist())
            elif len(nucleotides_above_lower) == 3:
                consensus += iupac_ambiguity(
                    nucleotides_above_lower.index.tolist())
            else:
                consensus += 'N'
        # Starts appending with the lower threshold
        else:
            if len(nucleotides_above_lower) == 3:
                consensus += iupac_ambiguity(
                    nucleotides_above_lower.index.tolist())
            elif len(nucleotides_above_middle) == 2:
                consensus += iupac_ambiguity(
                    nucleotides_above_middle.index.tolist())
            elif len(nucleotides_above_upper) == 1:
                consensus += nucleotides_above_upper.index[0]
            else:
                consensus += 'N'

    return consensus
