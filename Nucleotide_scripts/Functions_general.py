""" 
Script designed to recieve any length nucleotidic data, and create a file of PWM
"""

import sys
import os


def open_file(file_path: str) -> list:
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


def counter(nts_list: list) -> list[dict]:
    """Function that creates a list of dictionaries, that contains the count of nucleotides per position

    Args:
        nts_list (list): List of nucleotides

    Returns:
        list[dict]: Each item on the list represents the position, the dictionary contains the info of nucleotide count
    """

    nts_distribution = [{'A': 0, 'C': 0, 'G': 0, 'T': 0}
                        for n in range(len(nts_list[0]))]
    for line in nts_list:
        # for each "sequence" (line)
        for i in range(len(line)):
            # for every position, adds a count to the corresponding nts of that position
            nts_distribution[i][line[i]] += 1
    return nts_distribution


def pwmer(nts_list: list[dict], marker: str) -> str:
    """Function to return the Position Weight Matrix from the list of nucleotide count

    Args:
        nts_distribution (list[dict]): List of the count of nts per position
        marker (str): Marker to output in vertical or horizontal (-v/-h)

    Returns:
        str: A string of the PWM
    """

    if marker.lower() == "-h":
        nts_headers = ["A", "C", "G", "T"]
        final = f"PO    {'  '.join(map(str, range(1, len(nts_list) + 1)))}\n"
        for nts in range(4):
            # to traverse 4 times one for each nts
            nts_check = nts_headers.pop(0)
            line = f"{nts_check}    "
            for position in range(len(nts_list)):
                line += f"{str(nts_list[position][nts_check])}    "
            final += f"{line}\n"

    elif marker.lower() == "-v":
        final = "PO    A    C    G   T\n"
        for posicion in range(len(nts_list)):
            line = f"{str(posicion + 1)}    {nts_list[posicion]['A']}    {nts_list[posicion]['C']}    {nts_list[posicion]['G']}    {nts_list[posicion]['T']}"
            final += f"{line}\n"

    else:
        print("Please use a single letter V/H for vertical or horizontal ")
        return pwmer(nts_list)
    return final


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


def frecuency(nts_list: list) -> list[dict]:
    """Function to calculate the relative frequency of each nts per position

    Args:
        nts_list (list): List of nucleotides

    Returns:
        list[dict]: List where every item represents the position of the nucleotides and contains the info of relative frequency as a dictionary
    """
    nts_freq = [{'A': 0, 'C': 0, 'G': 0, 'T': 0}
                for n in range(len(nts_list))]
    for sequence in range(len(nts_list)):
        total = nts_list[sequence]["A"] + nts_list[sequence]["C"] + \
            nts_list[sequence]["G"] + nts_list[sequence]["T"]

        nts_freq[sequence]["A"] = round(nts_list[sequence]["A"] / total, 2)
        nts_freq[sequence]["C"] = round(nts_list[sequence]["C"] / total, 2)
        nts_freq[sequence]["G"] = round(nts_list[sequence]["G"] / total, 2)
        nts_freq[sequence]["T"] = round(nts_list[sequence]["T"] / total, 2)

    return nts_freq


def degen_iupac(nts: list) -> str:
    """Function to obtain the degenerate IUPAC nomenclature of the nts, according to the possible nts given

    Args:
        nts (list): List that contains only the possible possible nucleotides (in string format)

    Returns:
        str: A single letter of the corresponding degenerate IUPAC nomencalture
    """
    if len(nts) == 1:
        return nts[0]
    elif len(nts) == 2:
        if "A" in nts:
            if "T" in nts:
                return "W"
            elif "G" in nts:
                return "R"
            elif "C" in nts:
                return "M"
        elif "T" in nts:
            if "G" in nts:
                return "K"
            elif "C" in nts:
                return "Y"
        elif "G" in nts and "C" in nts:
            return "S"

    elif len(nts) == 3:
        if "A" in nts:
            if "T" in nts:
                if "G" in nts:
                    return "D"
                elif "C" in nts:
                    return "H"
            elif "C" in nts and "G" in nts:
                return "V"
        elif "G" in nts and "C" in nts and "T" in nts:
            return "B"


def ambigous_list(nts_freq: list[dict], threshold: str) -> list:
    """Function that uses the nts relative frequency and returns a list of the ambigous nts, per position

    Args:
        nts_freq (list[dict]): List of the relative frequency per position and nts
        threshold (str): Threshold to be considered when selecting posible nts

    Returns:
        list[list]: List containing lists of possible nucloteotides per position 
    """
    threshold = float(threshold)
    nts_list = []

    for pos in nts_freq:
        # Absolute nt
        # One item above the threshold and all other are bellow 0.2
        possible = filter(lambda nts: nts[1] > threshold and all(
            freq < 0.2 for key, freq in pos.items() if key != nts[0]), pos.items())
        nts = list(possible)
        if nts:
            nts_list.append([nts[0][0]])

        else:
            # 2 possible nts
            # if the condition above is not fulfilled
            # if one is above 62.5% and 2 of the rest are bellow 20%
            possible_above_threshold = filter(
                lambda nts: nts[1] > 0.625, pos.items())
            possible_inbetween_threshold = filter(
                lambda nts: 0.625 > nts[1] < 0.2, pos.items())
            above = list(possible_above_threshold)
            between = list(possible_inbetween_threshold)
            if len(above) == 1 and len(between) == 1:
                nts_list.append([above[0][0], between[0][0]])

            else:
                # 3 of possible nts
                # 1 is between 40 and 60%
                # 2the rest are above 20%
                possible_inbetween_threshold = filter(
                    lambda nts: 0.4 < nts[1] < 0.625, pos.items())
                possible_above_threshold = filter(
                    lambda nts: nts[1] > 0.2, pos.items())
                inbetween = list(possible_inbetween_threshold)
                above = list(possible_above_threshold)
                if len(above) == 2 and len(inbetween) == 1:
                    nts_list.append(
                        [above[0][0], above[1][0], inbetween[0][0]])
                else:
                    nts_list.append(['N'])
    return nts_list


def fuzzy_str(nts_list_ambigous) -> str:
    header = f"PO    {'  '.join(map(str, range(1, len(nts_list_ambigous) + 1)))}\n"
    nts_line = f"NTS"
    for pos in nts_list_ambigous:
        print(pos)
        nts = degen_iupac(pos)
        nts_line += f" {nts}"

    return header + nts_line


def regex_generator(sequence: list, max_mismatch: str):
    ### TO DO####
    # see the possibility of doing an alignment with regex and test if uses less pc
    regex_pattern = ""
    for pos in sequence:
        pos_pattern = ""
        for nucleotide_group in pos:
            if nucleotide_group == 'N':
                pos_pattern += '[ACGT]'
            else:
                pos_pattern += '[' + ''.join(nucleotide_group) + ']'
        regex_pattern += '(' + pos_pattern + ')' + \
            '{,' + str(max_mismatch) + '}'
    return regex_pattern


def alignment_brute(sequence: list, max_mismatch, target: str) -> list:
    """_summary_

    Args:
        sequence (list): _description_
        max_mismatch (_type_): _description_
        target (str): _description_

    Returns:
        _type_: _description_
    """
    succesfull_alignment = []
    max_mismatch = int(max_mismatch)
    for start in range(len(target)-len(sequence)):
        mismatch = 0
        excess = False
        print(start)
        for pos in range(len(sequence)):
            if sequence[pos][0] != 'N':
                if target[start+pos] not in sequence[pos]:

                    mismatch += 1
                    print(f"{target[start+pos]} {sequence[pos]}")
                    if mismatch > max_mismatch:
                        excess = True
                        break
        if excess == False:
            succesfull_alignment.append(start)

    return succesfull_alignment


if __name__ == "__main__":
    path_folder = "Input_files"
    try:
        file_name = sys.argv[1]
        analyze_marker = sys.argv[2]
        if analyze_marker == "-f":
            fuzzy_threshold = sys.argv[3]
            file_outputname = sys.argv[4]
        elif analyze_marker == "-p":
            orient_marker = sys.argv[3]
            file_outputname = sys.argv[4]
        elif analyze_marker == "-a":
            fuzzy_threshold = sys.argv[3]
            file_analyze = sys.argv[4]
            mismatch = sys.argv[5]
            file_outputname = sys.argv[6]
    except IndexError:
        print(
            "Please run the script with the name of the file to analyze, with the marker for PWM or fuzzy or both (-p/-f), \
if PWM, also with marker for horizontal or vertical (-h/-v), or if fuzzy, with threshold for recognition and with\
name(s) of the output file, please refer to the readme for more info"
        )

    try:
        path = os.path.join(path_folder, file_name)
        file = open_file(path)
        counted = counter(file)

    except FileNotFoundError:
        print("File name not found, please make sure to just add the name,and that the file to analyze is in the same folder \
as the script, refer to the readme for more info")
        sys.exit()
    if analyze_marker == "-p":
        if orient_marker == "-v" or orient_marker == "-h":
            output = pwmer(counted, orient_marker)
            freq = frecuency(counted)
            output2 = pwmer(freq, orient_marker)

    elif analyze_marker == "-f":
        freq = frecuency(counted)
        ambigous = ambigous_list(freq, fuzzy_threshold)
        output = fuzzy_str(ambigous)

    elif analyze_marker == "-a":
        freq = frecuency(counted)
        ambigous = ambigous_list(freq, fuzzy_threshold)
        analyze_path = os.path.join(path_folder, file_analyze)
        analyzed = open_file(analyze_path)
        analyzed = str(analyzed[0])
        aligned = alignment_brute(ambigous, mismatch, analyzed)
        output = str(aligned)
    else:
        print("Please use -p or -f for PWM or Fuzzy output")

    name = saver(output, f"{file_outputname}.abs")
    # if output2:
    #    saver(output2, f"{file_outputname}.rel")
    print(f"Succesfully created {name}")
