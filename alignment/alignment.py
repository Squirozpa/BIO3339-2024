"""
Functions to align query sequences to target sequences.
"""

# Standard Library Imports

# Local Library Imports

####################################################################################################


def is_match(query, target):
    """Internal function to check if a query nucleotide is a match for a target nucleotide."""
    # Aligment dictionary of degenerate nucleotides (IUPAC codes)
    alignment_dict = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG",
                      "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}

    # Get the possible matches for the query
    possible_matches = alignment_dict.get(query, "")
    # Check if the target is a possible match
    return any(char in possible_matches for char in target)


def reverse_complement(sequence: str) -> str:
    """Internal function te generate the reverse complement of a DNA sequence."""
    swap_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join([swap_dict[char] for char in sequence[::-1]])


def align(query_sequence: str, target_sequence: str, max_mismatch: int = 0) -> dict[int, list]:
    """
    Aligns a query sequence with a target sequence and returns a dictionary of alignments.

    Args:
        query_sequence (str): The query sequence to align.
        target_sequence (str): The target sequence to align against.
        max_mismatch (int, optional): The maximum number of allowed mismatches. Defaults to 0.

    Returns:
        dict[int, list]: A dictionary where the keys represent the number of mismatches and the values are lists of start positions.

    """
    # Aligment dictionary of degenerate nucleotides (IUPAC codes)
    alignment_dict = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG",
                      "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
    # Dictionary to store the alignments
    alignments = {}
    # Starts iterating over the length of the target sequence minus the length of the query sequence
    for start in range(len(target_sequence) - len(query_sequence) + 1):
        # Starts with 0 mismatches
        mismatches = 0
        # Iterates over the length of the query sequence
        for position in range(len(query_sequence)):
            # If the query sequence nucleotide is not in the alignment dictionary of the target sequence nucleotide
            if not is_match(query_sequence[position], [target_sequence[start + position]]):
                # Adds 1 to the mismatches
                mismatches += 1
            # If the mismatches are greater than the maximum allowed mismatches, breaks the loop
            if mismatches > max_mismatch:
                break
        if mismatches <= max_mismatch:
            if mismatches in alignments:
                # If theres is a key for that mismatch appends to the key of the mismatches, the start position
                alignments[mismatches].append(start)
            else:
                # If there is no key for that mismatch, creates a key for that mismatch and creates a list with the start position
                alignments[mismatches] = [start]
    return alignments


def align_both_directions(query_sequence: str, target_sequence: str, max_mismatch: int = 1) -> list[list]:
    """
    Aligns the query sequence in both the original and reversed direction of the target sequence.

    Args:
        query_sequence (str): The sequence to be aligned.
        target_sequence (str): The target sequence to align against.
        max_mismatch (int, optional): The maximum number of allowed mismatches. Defaults to 1.

    Returns:
        list[list]: A list containing two sublists. The first sublist contains the matches found in the original target sequence,
                    and the second sublist contains the matches found in the reversed target sequence.
    """
    # Find matches in the original target sequence
    matches_original = align(query_sequence, target_sequence, max_mismatch)

    # Find the reverse complement of the target sequence
    reversed_sequence = reverse_complement(target_sequence)
    # Find matches in the reversed target sequence
    matches_reversed = align(query_sequence, reversed_sequence, max_mismatch)
    # Change all values in the dictionary
    for key in matches_reversed:
        matches_reversed[key] = [len(
            target_sequence) - value - 1 for value in matches_reversed[key]]
    # Combine the results
    return [matches_original, matches_reversed]
