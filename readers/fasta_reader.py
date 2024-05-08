"""
Functions to create a string for exporting a matrix to a tab-separated file
"""

# Standard Library Imports

# Local Library Imports

####################################################################################################

"""
Function to save the sequences in a fasta file as a list of strings.
"""

# Standard Library Imports

# Local Library Imports

####################################################################################################


def read_fasta(fasta_file_path: str, raw: bool = False) -> list[str]:
    """
    Read a fasta file and return each sequence as a str in a list.
    Args:
        fasta_file_path (str): Path to the fasta file.
        raw (bool, optional): If the file is a .seq file. Defaults to False.
    Returns:
        list[str]: List of nucleotide sequences in the fasta file.
    """
    sequence_list = []
    sequence = ""
    with open(fasta_file_path, "r") as fasta_file:
        if raw:
            for line in fasta_file:
                sequence_list.append(line.strip())
        else:
            for line in fasta_file:
                if line.startswith(">"):
                    if sequence:
                        sequence_list.append(sequence)
                        sequence = ""
                else:
                    sequence += line.strip()
            if sequence:
                sequence_list.append(sequence)  # Save the last sequence
    return sequence_list


if __name__ == "__main__":
    print("This file is not meant to be run directly.")
