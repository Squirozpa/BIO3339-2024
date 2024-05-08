"""
Functions to export the alignment to a file
"""

# Standard Library Imports

# Local Library Imports

####################################################################################################


def write_alignment(alignment: list[dict], output_file: str, max_mismatch: int = None) -> str:
    """
    Function to write the alignment to a file and save it as a txt file.

    Parameters:
    - alignment (list[dict]): The alignment as a list of dictionaries.
    - output_file (str): The path to the output file.

    Returns:
    - output_file (str): The path to the output file.

    """
    titles = ["Alignment Wattson", "Alignment Crick"]
    with open(f"{output_file}.txt", 'w') as f:
        for i, align_dict in enumerate(alignment):
            # Saves the Wattson alignment (first item on the list) then Crick alignment (second item on the list)
            f.write(f"{titles[i]}:\n")
            for key, value in align_dict.items():
                # Saves the mismatches for each alignment
                f.write(f"Mismatches({key}): {value}\n")
            f.write("\n")
        if max_mismatch:
            # Adds the max mismatches to the file
            f.write(f"Max Mismatches: {str(max_mismatch)}")
    return output_file
