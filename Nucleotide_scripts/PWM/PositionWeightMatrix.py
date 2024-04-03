"""Collection of functions associated with the creation of the position weight matrix"""
import math


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


def pwmer(nts_list: list[dict]) -> str:
    """Function to return the Position Weight Matrix from the list of nucleotide count

    Args:
        nts_distribution (list[dict]): List of the count of nts per position
        marker (str): Marker to output in vertical or horizontal (-v/-h) ###marker use is discontinued

    Returns:
        str: A string of the PWM
    """

    nts_headers = ["A", "C", "G", "T"]
    final = f"PO    {'  '.join(map(str, range(1, len(nts_list) + 1)))}\n"
    for nts in range(4):
        # to traverse 4 times one for each nts
        nts_check = nts_headers.pop(0)
        line = f"{nts_check}    "
        for position in range(len(nts_list)):
            line += f"{str(nts_list[position][nts_check])}    "
        final += f"{line}\n"
    """
    elif marker.lower() == "-v":
        final = "PO    A    C    G   T\n"
        for posicion in range(len(nts_list)):
            line = f"{str(posicion + 1)}    {nts_list[posicion]['A']}    {nts_list[posicion]['C']}    {nts_list[posicion]['G']}    {nts_list[posicion]['T']}"
            final += f"{line}\n"
    """
    """else:
        print("Please use a single letter V/H for vertical or horizontal ")
        return pwmer(nts_list)
    """
    return final


def frequency(nts_list: list) -> list[dict]:
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


def weight(nts_freq: list[dict]) -> list[dict]:
    """Function to calculate the weight of each nts per position

    Args:
        nts_freq (list[dict]): List of the relative frequency per position and nts

    Returns:
        list[dict]: List where every item represents the position of the nucleotides and contains the info of weight as a dictionary
    """
    nts_weight = [{'A': 0, 'C': 0, 'G': 0, 'T': 0}
                  for n in range(len(nts_freq))]
    inf_values = 0.0000000001
    for sequence in range(len(nts_freq)):
        total = nts_freq[sequence]["A"] + nts_freq[sequence]["C"] + \
            nts_freq[sequence]["G"] + nts_freq[sequence]["T"]
        if nts_freq[sequence]["A"] == 0:
            nts_freq[sequence]["A"] = inf_values
        if nts_freq[sequence]["C"] == 0:
            nts_freq[sequence]["C"] = inf_values
        if nts_freq[sequence]["G"] == 0:
            nts_freq[sequence]["G"] = inf_values
        if nts_freq[sequence]["T"] == 0:
            nts_freq[sequence]["T"] = inf_values

        nts_weight[sequence]["A"] = round(math.log(
            nts_freq[sequence]["A"] / 0.25, 2), 2)
        nts_weight[sequence]["C"] = round(math.log(
            nts_freq[sequence]["C"] / 0.25, 2), 2)
        nts_weight[sequence]["G"] = round(math.log(
            nts_freq[sequence]["G"] / 0.25, 2), 2)
        nts_weight[sequence]["T"] = round(math.log(
            nts_freq[sequence]["T"] / 0.25, 2), 2)
    return nts_weight


def max_score_prm(PRM: list[dict]) -> float:
    """Function to calculate the max score of a sequence, given a Position Weight Matrix

    Args:
        PXM (list[dict]): List of the weight of each nts per position

    Returns:
        float: Max core of the sequence
    """
    max_score = 1
    for item in PRM:
        max_score *= max(item.values())
    return max_score


def score_prm(PRM: list[dict], sequence: str) -> list[list, float]:
    """Function to calculate the score of a sequence, given a Position Weight Matrix

    Args:
        PXM (list[dict]): List of the weight of each nts per position
        sequence (str): Sequence to calculate the score

    Returns:
        float: Score of the sequence
    """
    score_list = []
    score = 1
    for pos in range(len(PRM)):
        score *= PRM[pos][sequence[pos]]
        score_list.append(score)
    return [score_list, score]


def max_score_pwm(PWM: list[dict]) -> float:
    """Function to calculate the max score of a sequence, given a Position Weight Matrix

    Args:
        PXM (list[dict]): List of the weight of each nts per position

    Returns:
        float: Max core of the sequence
    """
    max_score = 0
    for item in PWM:
        max_score += max(item.values())
    return max_score


def score_pwm(PWM: list[dict], sequence: str) -> list[list, float]:
    """Function to calculate the score of a sequence, given a Position Weight Matrix

    Args:
        PXM (list[dict]): List of the weight of each nts per position
        sequence (str): Sequence to calculate the score

    Returns:
        float: Score of the sequence
    """
    score_list = []
    score = 0
    for pos in range(len(PWM)):
        score += PWM[pos][sequence[pos]]
        score_list.append(score)
    return [score_list, score]
