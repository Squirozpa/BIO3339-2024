"""
Function to export a sequence in FASTA format
"""

# Standard Library Imports
import os
# Local Library Imports

####################################################################################################


def write_fasta(sequence: str, sequence_id: str, output_file: str) -> str:
    """
    Function to export a sequence in FASTA format and save it as a .fasta file.

    Args:
        sequence (str): Sequence to be exported
        sequence_id (str): ID of the sequence
        output_file (str): Name of the output file

    Returns:
        output_file (str): Name of the output file
    """
    # Create the FASTA string
    # Limiting the length of the sequence to 80 characters
    sequence_lines = [sequence[i:i+80] for i in range(0, len(sequence), 80)]
    fasta_str = ">" + sequence_id+"\n" + "\n".join(sequence_lines) + "\n"

    # Save the FASTA string as a .fasta file
    with open(f"{output_file}.fasta", 'w') as file:
        file.write(fasta_str)

    return output_file
