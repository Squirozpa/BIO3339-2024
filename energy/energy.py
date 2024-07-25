"""
Functions to calculate the energy of an sequence in a sequence, from a matrix
"""

# Standard Library Imports
import pandas as pd
# Local Library Imports
from matrices import scores as ms
####################################################################################################


def reverse_complement(sequence: str) -> str:
    """Internal function te generate the reverse complement of a DNA sequence."""
    swap_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join([swap_dict[char] for char in sequence[::-1]])


def genomic_energy(matrix: pd.DataFrame, sequence: str, matrix_type: str = "weight", threshold: float = None, normalized: str = "normalized", reversed: bool = True, only_reversed: bool = False) -> pd.DataFrame:
    """
    Function that calculates the energy of a sequence in a genomic sequence, given a matrix.

    Parameters:
    - matrix (pd.DataFrame): The weight or relative matrix used for energy calculation.
    - sequence (str): The genomic sequence for which energy needs to be calculated.
    - matrix_type (str, optional): The type of matrix used for energy calculation. Default is "weight".
    - threshold (float, optional): Minimum energy value. If an energy value is lower than this, it will be set to this value. Default is None.
    - normalized (str, optional): Flag indicating whether the energy values should be normalized. Options are "normalized", "relative", or "none". Default is "normalized".
    - reversed (bool, optional): Flag indicating whether the reversed energy values should also be calculated. Default is True.
    - only_reversed (bool, optional): Flag indicating whether only the reversed energy values should be calculated. Default is False.

    Returns:
    - pd.DataFrame: A DataFrame of energy values for the given sequence. If `reversed` is True, it also returns reversed energy values.

    Raises:
    - ValueError: If the `matrix_type` is not "weight" or "relative".

    """
    print(reversed, only_reversed)
    if reversed == False and only_reversed:
        raise ValueError(
            "--no_reversed and --only_reversed cannot be used together")
    # Initialize the DataFrame of energies
    genomic_energy = pd.DataFrame(
        columns=['energy']) if not only_reversed else None
    # If reversed is needed, also initialize the DataFrame of reversed energies
    if reversed or only_reversed:
        reversed_sequence = reverse_complement(sequence)
        reversed_genomic_energy = pd.DataFrame(columns=['reversed_energy'])
    ##### Block for weight matrix, used to reduce if statements when running#####
    if matrix_type == "weight":
        # Calculate the max and min scores for normalization if needed
        if normalized == "relative":
            max_score = ms.max_score_weight_matrix(matrix)
            min_score = ms.min_score_weight_matrix(matrix)
        elif normalized == "normalized":
            max_score = ms.max_score_weight_matrix(matrix)
        elif normalized == "none":
            pass
        else:
            raise ValueError(
                "Normalization must be 'normalized', 'relative', or 'none'")
        # Iterate over the sequence and calculate the energy
        for i in range(len(sequence) - len(matrix) + 1):
            # Calculates the energy of the sequence in that position
            if not only_reversed:
                energy = ms.score_weight_matrix(
                    sequence[i:i + len(matrix)], matrix)
            # Also calculates the reversed energy if needed
            if reversed:
                reverse_energy = ms.score_weight_matrix(
                    reversed_sequence[i:i + len(matrix)], matrix)
                # If threshold is provided and energy is lower than threshold, set energy to threshold
                if threshold is not None and reverse_energy < threshold:
                    reverse_energy = threshold
                # If normalized, divide by the max score
                elif normalized == "relative":
                    reverse_energy = ms.relative_score(
                        min_score, max_score, reverse_energy)
                elif normalized == "normalized":
                    reverse_energy = reverse_energy / max_score
                elif normalized == "none":
                    pass
            # If threshold is provided and energy is lower than threshold, set energy to threshold
            if not only_reversed:
                if threshold is not None and energy < threshold:
                    energy = threshold
                # If normalized, divide by the max score
                elif normalized == "relative":
                    energy = ms.relative_score(min_score, max_score, energy)
                elif normalized == "normalized":
                    energy = energy / max_score
                elif normalized == "none":
                    pass
            # Append the energy to the DataFrame if not only_reversed
                genomic_energy = genomic_energy._append(
                    {'energy': energy}, ignore_index=True)
            if reversed:
                # Also append the reversed energy to the DataFrame if needed
                reversed_genomic_energy = reversed_genomic_energy._append(
                    {'reversed_energy': reverse_energy}, ignore_index=True)
        if reversed:
            # Concatenate the two DataFrames if reversed energies are needed
            reversed_genomic_energy = reversed_genomic_energy.iloc[::-1].reset_index(
                drop=True)
            genomic_energy = pd.concat(
                [genomic_energy, reversed_genomic_energy], axis=1)
        if not only_reversed and reversed:
            genomic_energy = reversed_genomic_energy.iloc[::-1].reset_index(
                drop=True)
    ##### Block for relative matrix, used to reduce if statements when running#####
    elif matrix_type == "relative":
        # Calculate the max and min scores for normalization if needed
        if normalized == "relative":
            max_score = ms.max_score_weight_matrix(matrix)
            min_score = ms.min_score_weight_matrix(matrix)
        elif normalized == "normalized":
            max_score = ms.max_score_weight_matrix(matrix)
        elif normalized == "none":
            pass
        else:
            raise ValueError(
                "Normalization must be 'normalized', 'relative', or 'none'")
        # Iterate over the sequence and calculate the energy
        for i in range(len(sequence) - len(matrix) + 1):
            # Calculates the energy of the sequence in that position
            energy = ms.score_weight_matrix(sequence[i:i + len(matrix)], matrix)
            # Also calculates the reversed energy if needed
            if reversed:
                reverse_energy = ms.score_weight_matrix(
                    reversed_sequence[i:i + len(matrix)], matrix)
                # If threshold is provided and energy is lower than threshold, set energy to threshold
                if threshold is not None and reverse_energy < threshold:
                    reverse_energy = threshold
                # If normalized, divide by the max score
                elif normalized == "relative":
                    reverse_energy = ms.relative_score(
                        min_score, max_score, reverse_energy)
                elif normalized == "normalized":
                    reverse_energy = reverse_energy / max_score
                elif normalized == "none":
                    pass
            # If threshold is provided and energy is lower than threshold, set energy to threshold
            if threshold is not None and energy < threshold:
                energy = threshold
            # If normalized, divide by the max score
            elif normalized == "relative":
                energy = ms.relative_score(min_score, max_score, energy)
            elif normalized == "normalized":
                energy = energy / max_score
            elif normalized == "none":
                pass
            # Append the energy to the DataFrame
            genomic_energy = genomic_energy._append(
                {'energy': energy}, ignore_index=True)
            if reversed:
                # Also append the reversed energy to the DataFrame if needed
                reversed_genomic_energy = reversed_genomic_energy._append(
                    {'reversed_energy': reverse_energy}, ignore_index=True)
        if reversed:
            # Concatenate the two DataFrames if reversed energies are needed
            reversed_genomic_energy = reversed_genomic_energy.iloc[::-1].reset_index(
                drop=True)
            genomic_energy = pd.concat(
                [genomic_energy, reversed_genomic_energy], axis=1)
    else:
        raise ValueError("Matrix type must be 'weight' or 'relative'")

    return genomic_energy
