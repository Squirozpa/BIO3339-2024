"""
Function to export a matrix to a file
"""

# Standard Library Imports
import pandas as pd
# Local Library Imports

####################################################################################################


def write_matrix(matrix: pd.DataFrame, output_file: str, max_value: float = None, vertical: bool = False, decimals: bool = False) -> str:
    """
    Function to export a matrix in tab-separated format and save it as a .tsv file.

    Args:
        matrix (pd.DataFrame): Matrix to be exported
        output_file (str): Path to save the output file
        max_value (float, optional): Maximum value to be added at the end of the file. Defaults to None.
        vertical (bool, optional): Flag to indicate if the matrix should be exported vertically. Defaults to False.
        decimals (bool, optional): Flag to indicate if the matrix should be exported with 2 decimals. Defaults to False.
    Returns:
        str: Name of the output file
    """
    if decimals:
        # Convert the DataFrame to a string with tab separation
        matrix = matrix.round(3)  # Round to 2 decimal places
        # Format as string with 2 decimal places
        matrix = matrix.map('{:.2f}'.format)

    if vertical:
        matrix_str = matrix.to_csv(
            sep="\t", index=True, header=True, index_label="PO")
    else:
        matrix_str = matrix.T.to_csv(
            sep="\t", index=True, header=True, index_label="PO")

    # Add max value if present
    if max_value is not None:
        matrix_str += f"#max value: {max_value}"

    # Save the matrix to the output file
    with open(f"{output_file}.tsv", "w") as file:
        file.write(matrix_str)

    return output_file
