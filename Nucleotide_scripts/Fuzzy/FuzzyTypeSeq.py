"""Collection of functions associated with the creation of the power weight matrix"""


def degen_iupac(nts_unique: list) -> str:
    """Function to obtain the degenerate IUPAC nomenclature of the nts, according to the possible nts given

    Args:
        nts (list): List that contains only the possible possible nucleotides (in string format)

    Returns:
        str: A single letter of the corresponding degenerate IUPAC nomencalture
    """
    if len(nts_unique) == 1:
        return nts_unique[0]
    elif len(nts_unique) == 2:
        if "A" in nts_unique:
            if "T" in nts_unique:
                return "W"
            elif "G" in nts_unique:
                return "R"
            elif "C" in nts_unique:
                return "M"
        elif "T" in nts_unique:
            if "G" in nts_unique:
                return "K"
            elif "C" in nts_unique:
                return "Y"
        elif "G" in nts_unique and "C" in nts_unique:
            return "S"

    elif len(nts_unique) == 3:
        if "A" in nts_unique:
            if "T" in nts_unique:
                if "G" in nts_unique:
                    return "D"
                elif "C" in nts_unique:
                    return "H"
            elif "C" in nts_unique and "G" in nts_unique:
                return "V"
        elif "G" in nts_unique and "C" in nts_unique and "T" in nts_unique:
            return "B"
    else:
        return "N"


def ambigous_nts(nts_freq: list[dict], threshold: str) -> list:
    """Creates a list of the ambigous nts, per position

    Args:
        nts_freq (dict): List of the relative frequency per position and nts
        threshold (str): Threshold to be considered when selecting posible nts

    Returns:
        list: List containing the possible nucloteotides per position
    """
    threshold = float(threshold)
    list_ambigous_nts = []
    for pos in range(len(nts_freq)):
        max_value = float(max(nts_freq[pos].values()))

        effective_treshold = threshold*max_value
        filtered = filter(
            lambda nts: float(nts[1]) >= effective_treshold, nts_freq[pos].items())
        to_append = dict(filtered)
        list_ambigous_nts.append(list(to_append.keys()))

    return list_ambigous_nts


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


def fuzzy_str(nts_list_ambigous, threshold) -> str:
    """Converts the list the the string type for the fuzzy type sequence

    Args:
        nts_list_ambigous (_type_): List of ambigous nucleotides in following IUPAC convention

    Returns:
        str: A 2 line string for saving the Fuzzy type sequence
    """
    header = "PO\t" + \
        "\t".join(map(str, range(1, len(nts_list_ambigous) + 1))) + "\n"
    nts_line = f"NTS"
    for pos in nts_list_ambigous:
        nts = degen_iupac(pos)
        nts_line += "\t" + nts
    nts_line += "\n"+f"#threshold: {threshold}"

    return header + nts_line
