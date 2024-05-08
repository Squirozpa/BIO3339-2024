"""
Function to save matrix as a pandas DataFrame. Every row represents the position of the sequence
"""

# Standard Library Imports
import pandas as pd
# Local Library Imports

####################################################################################################


def read_matrix(input_file: str) -> pd.DataFrame:
    """
    Function to read a matrix from a .tsv file and convert it into a DataFrame.

    Args:
        input_file (str): Path to the input file

    Returns:
        pd.DataFrame: DataFrame with the matrix
    """
    # Read the file and convert it into a DataFrame
    matrix = pd.read_csv(input_file, sep="\t", index_col="PO", comment='#')

    # If the matrix is horizontal, transpose it to make it vertical
    if matrix.shape[0] < matrix.shape[1]:
        matrix = matrix.T

    # Convert the data to numbers
    matrix = matrix.apply(pd.to_numeric, errors='coerce')

    # Reset the index and drop the 'PO' column
    matrix.reset_index(drop=True, inplace=True)

    return matrix
