"""
Script designed to recieve any length nucleotidic data, and create a file of PWM
"""

import sys
import os
from collections import defaultdict

import Fuzzy.FuzzyTypeSeq as fuzzy
import PWM.PositionWeightMatrix as pwm
import Alignment.Alignment as ali


def open_fasta(file_path: str) -> list:
    """Function that recives the path to open the file then returns a list of the nts

    Args:
        file_path (str): File name for processing

    Returns:
        list: List of nts, where every line is a different sequence
    """

    file = open(file_path, "r")
    lines_file = file.read().splitlines(False)
    file.close()
    nts_list = []
    for line in lines_file:
        if line.startswith(">"):
            pass
        else:
            nts_list.append(line)
    return nts_list


def open_pwm(file_path: str) -> list[dict]:
    """Function that recives the path to open the file then returns a list of the value per nts of a pwm type table

    Args:
        file_path (str): File name for processing

    Returns:
        list: list of dictionaries, where every dictionary is a position of the nts
    """

    file = open(file_path, "r")
    lines_file = file.read().splitlines(False)
    file.close()
    file_lines = []
    for line in lines_file:
        if line.startswith("PO"):
            pass
        else:
            line = line.split("\t")[1:]
            file_lines.append(line)
    freq_list = [{'A': 0, 'C': 0, 'G': 0, 'T': 0}
                 for n in range(len(file_lines[0]))]
    for pos in range(len(file_lines[0])):
        freq_list[pos]["A"] = float(file_lines[0][pos])
        freq_list[pos]["C"] = float(file_lines[1][pos])
        freq_list[pos]["G"] = float(file_lines[2][pos])
        freq_list[pos]["T"] = float(file_lines[3][pos])

    return freq_list


def open_fuzzy(file_path: str) -> list:
    """Function that recives the path to open the file then returns a list of the value per nts of a fuzzy type table
    Args:
        file_path (str): File name for processing
    Returns:
        list: list of ambigous nts
    """

    file = open(file_path, "r")
    lines_file = file.read().splitlines(False)
    file.close()
    fuzzy_list = lines_file[1].split("\t")[1:]
    return fuzzy_list


def saver(save_string: str, file_name: str) -> str:
    """Function to save the file as a .txt

    Args:
        save_string (str): String containing the info to save as a file
        file_name (str): Name of the file to be saved as

    Returns:
        str: File name
    """

    file = open(f"Output_files/{file_name}", "w")
    file.write(save_string)
    return file_name


def pwm_generator(input_path: str, ouput_name: str) -> None:
    """Function to generate the PWM from the input file

    Args:
        raw_input_path (str): Name of the file containing the raw data in the input folder
        ouput_name (str): Name of the file to save the output tables
    """
    if input_path.startswith("="):
        raw_input_path = os.path.join(
            os.getcwd(), "Output_files", input_path[1:])
    else:
        raw_input_path = os.path.join(os.getcwd(), "Input_files", input_path)
    nts_list = open_fasta(raw_input_path)
    list_count = pwm.pfm(nts_list)
    table_count = pwm.matrix_writer(list_count, 0)
    list_freq = pwm.ppm(list_count)
    max_score_ppm = pwm.max_score_ppm(list_freq)
    table_freq = pwm.matrix_writer(list_freq, max_score_ppm)
    list_weight = pwm.pwm(list_freq)
    max_score_pwm = pwm.max_score_pwm(list_weight)
    table_weight = pwm.matrix_writer(list_weight, max_score_pwm)
    saver(table_count, f"{ouput_name}_pfm.tab")
    saver(table_freq, f"{ouput_name}_ppm.tab")
    saver(table_weight, f"{ouput_name}_pwm.tab")

    return None


def fuzzy_generator(input_path: str, ouput_name: str, treshold: str) -> None:
    """
    Generates a fuzzy string based on the input file.

    Parameters:
    input_path (str): The path to the input file.
    ouput_name (str): The name of the output file.
    treshold (str): The threshold value for determining ambiguous nucleotides.

    Returns:
    None
    """
    if input_path.startswith("="):
        freq_input_path = os.path.join(
            os.getcwd(), "Output_files", input_path[1:])
    else:
        freq_input_path = os.path.join(os.getcwd(), "Input_files", input_path)
    nts_freq = open_pwm(freq_input_path)
    ambigous_list = fuzzy.ambigous_nts(nts_freq, treshold)
    # ambigous_list = fuzzy.ambigous_list(nts_freq, treshold)
    ambigous_str = fuzzy.fuzzy_str(ambigous_list, treshold)
    saver(ambigous_str, f"{ouput_name}_fuzzy.tab")

    return None


def save_alignment(alignment, output_name):
    sting_to_save = "Alignment Wattson:\n"
    grouped_alignments = defaultdict(list)
    for position, mismatches in alignment[0].items():
        grouped_alignments[mismatches].append(position)
    for mismatches in sorted(grouped_alignments.keys()):
        positions = ', '.join(map(str, sorted(grouped_alignments[mismatches])))
        sting_to_save += f"Mismatches: {mismatches}, Positions: {positions}\n"
    sting_to_save += "Alignment Crick:\n"
    grouped_alignments = defaultdict(list)
    for position, mismatches in alignment[1].items():
        grouped_alignments[mismatches].append(position)
    for mismatches in sorted(grouped_alignments.keys()):
        positions = ', '.join(map(str, sorted(grouped_alignments[mismatches])))
        sting_to_save += f"Mismatches: {mismatches}, Positions: {positions}\n"
    with open(f"Output_files/{output_name}_ali.tab", "w") as file:
        file.write(sting_to_save)
    return None


def alignment_generator(input_path_seq: str, input_path_target: str, ouput_name: str, max_mismatch: str) -> None:
    """
    Generates an alignment using the given input sequence and target file.

    Args:
        input_path_seq (str): The name to the input sequence file.
        input_path_target (str): The name to the target file.
        ouput_name (str): The name of the output file.
        max_mismatch (str): The maximum number of allowed mismatches.

    Returns:
        None
    """
    if input_path_seq.startswith("="):
        fuzzy_input_path = os.path.join(
            os.getcwd(), "Output_files", input_path_seq[1:])
    else:
        fuzzy_input_path = os.path.join(
            os.getcwd(), "Input_files", input_path_seq)
    fuzzy_list = open_fuzzy(fuzzy_input_path)
    if input_path_target.startswith("="):
        input_path_target = os.path.join(
            os.getcwd(), "Output_files", input_path_target[1:])
    else:
        input_path_target = os.path.join(
            os.getcwd(), "Input_files", input_path_target)
    target = open_fasta(input_path_target)
    max_mismatch = int(max_mismatch)
    alignment = ali.alignment_brute(
        target=target[0], sequence=fuzzy_list, max_mismatch=max_mismatch)
    save_alignment(alignment, ouput_name)

    return None


def energy_writer(energy_values: list) -> str:
    """
    Writes the energy values to a string.

    Args:
        energy_values (list): The list of energy values.

    Returns:
        str: The string containing the energy values.
    """
    energy_string = '\t'.join(str(value) for value in energy_values[0])+"\n"
    energy_string += '\t'.join(str(value) for value in energy_values[1])+"\n"

    return energy_string


def energy_generator(prm_pwm_path: str, genome_path: str, prm_pwm_type: str, output_name: str) -> None:
    """
    Generates energy values based on the given prm/pwm and genome.

    Args:
        prm_pwm_path (str): The path to the prm/pwm file.
        genome_path (str): The path to the genome file.
        prm_pwm_type (str): The type of prm/pwm ('prm' or 'pwm').
        create_plot (bool, optional): Whether to create a plot or not. Defaults to False.

    Returns:
        None
    """
    # Read prm/pwm file
    if prm_pwm_path.startswith("="):
        prm_pwm_path = os.path.join(
            os.getcwd(), "Output_files", prm_pwm_path[1:])
    else:
        prm_pwm_path = os.path.join(os.getcwd(), "Input_files", prm_pwm_path)
    prm_pwm_list = open_pwm(prm_pwm_path)

    # Read genome file
    if genome_path.startswith("="):
        genome_path = os.path.join(
            os.getcwd(), "Output_files", genome_path[1:])
    else:
        genome_path = os.path.join(os.getcwd(), "Input_files", genome_path)
    genome_list = open_fasta(genome_path)
    # get the energy values
    values = ali.genomic_energy_profile(
        matrix=prm_pwm_list, genome=genome_list[0], matrix_type=prm_pwm_type)

    # Save energy values to file
    energy_string = energy_writer(values)
    saver(energy_string, f"{output_name}_energy.tab")

    return None


if __name__ == "__main__":
    try:
        if sys.argv[1] == "pwm":
            pwm_generator(sys.argv[2], sys.argv[3])
        elif sys.argv[1] == "fuzzy":
            fuzzy_generator(sys.argv[2], sys.argv[3], sys.argv[4])
        elif sys.argv[1] == "alignment":
            alignment_generator(
                sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
        elif sys.argv[1] == "energy":
            energy_generator(sys.argv[2], sys.argv[3],
                             sys.argv[4], sys.argv[5])
        else:
            print("Invalid function")
    except IndexError:
        print("Invalid arguments, please refer to the README for more information")
    sys.exit(1)
