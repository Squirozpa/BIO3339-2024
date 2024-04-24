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
        nts_in_pos = []

        # Only 1 nts
        for nts in nts_freq[pos]:
            if nts_freq[pos][nts] >= (0.25 + threshold*(0.75)):
                nts_in_pos.append(nts)

        if not nts_in_pos:
            # 2 nts
            for nts in nts_freq[pos]:
                if nts_freq[pos][nts] >= (0.25 + threshold*(0.25)):
                    nts_in_pos.append(nts)
            if len(nts_in_pos) != 2:
                nts_in_pos = []

            # 3 nts
        if not nts_in_pos:
            for nts in nts_freq[pos]:
                if nts_freq[pos][nts] >= (0.25 + threshold*(0.08)):
                    nts_in_pos.append(nts)

            if len(nts_in_pos) != 3:
                nts_in_pos = []
        print(nts_in_pos)
        if nts_in_pos:
            list_ambigous_nts.append(degen_iupac(nts_in_pos))
        else:
            list_ambigous_nts.append("N")
    return list_ambigous_nts


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
