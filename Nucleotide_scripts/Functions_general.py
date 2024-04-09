""" 
Script designed to recieve any length nucleotidic data, and create a file of PWM
"""

import sys
import os

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
            line = line.split()[1:]
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
    fuzzy_list = lines_file[1].split()[1:]
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
    raw_input_path = os.path.join(os.getcwd(), "Input_files", input_path)
    nts_list = open_fasta(raw_input_path)
    list_count = pwm.pfm(nts_list)
    table_count = pwm.matrix_writer(list_count, 0)
    list_freq = pwm.ppm(list_count)
    max_score_prm = pwm.max_score_pwm(list_freq)
    table_freq = pwm.matrix_writer(list_freq, max_score_prm)
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

    freq_input_path = os.path.join(os.getcwd(), "Input_files", input_path)
    nts_freq = open_pwm(freq_input_path)
    ambigous_list = fuzzy.ambigous_nts(nts_freq, treshold)
    # ambigous_list = fuzzy.ambigous_list(nts_freq, treshold)
    ambigous_str = fuzzy.fuzzy_str(ambigous_list)
    saver(ambigous_str, f"{ouput_name}.fuzzy")

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
    fuzzy_input_path = os.path.join(os.getcwd(), "Input_files", input_path_seq)
    fuzzy_list = open_fuzzy(fuzzy_input_path)
    input_path_target = os.path.join(
        os.getcwd(), "Input_files", input_path_target)
    target = open_fasta(input_path_target)
    max_mismatch = int(max_mismatch)
    alignment = ali.alignment_brute(
        target=target[0], sequence=fuzzy_list, max_mismatch=max_mismatch)
    saver(str(alignment), f"{ouput_name}.ali")


if __name__ == "__main__":
    if sys.argv[1] == "pwm":
        pwm_generator(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "fuzzy":
        fuzzy_generator(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "alignment":
        alignment_generator(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    else:
        print("Invalid function")
