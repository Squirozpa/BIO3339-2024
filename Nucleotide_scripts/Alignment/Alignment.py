"""Collection of functions associated with the alignment    
"""


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
