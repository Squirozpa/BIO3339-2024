"""
Funtion to export genomic energy profile to a file
"""

# Standard Library Imports
import pandas as pd
# Local Library Imports

####################################################################################################


def write_energy(genomic_energy: pd.DataFrame, output_file: str) -> str:
    """
    Function to export the genomic profile to a file and save it as a .tsv file.

    Parameters:
    - genomic_energy (pd.DataFrame): The DataFrame of energy values for the genomic sequence.
    - output_file (str): The path to the output file.

    Returns:
    - output_file (str): The path to the output file.

    """
    genomic_energy.to_csv(output_file, sep='\t', index=False)

    return output_file
