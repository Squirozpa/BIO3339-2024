"""Collection of functions associated with the alignment
"""
import sys
import PWM.PositionWeightMatrix as pwm
import matplotlib.pyplot as plt


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
    alignment_dict = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG",
                      "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
    succesfull_alignment_wattson = {}
    succesfull_alignment_crick = {}
    max_mismatch = int(max_mismatch)
    for start in range(len(target)-len(sequence)+1):
        mismatch = 0
        excess = False
        for pos in range(len(sequence)):
            if sequence[pos][0] != 'N':
                if target[start+pos] not in alignment_dict[sequence[pos]]:
                    mismatch += 1
                    if mismatch > max_mismatch:
                        excess = True
                        break
        if excess == False:
            succesfull_alignment_wattson[start+1] = mismatch
    reversed_target = target[::-1]
    for start in range(len(target)-len(sequence)+1):
        mismatch = 0
        excess = False
        for pos in range(len(sequence)):
            if sequence[pos][0] != 'N':
                if reversed_target[start+pos] not in alignment_dict[sequence[pos]]:
                    mismatch += 1
                    if mismatch > max_mismatch:
                        excess = True
                        break
        if excess == False:
            succesfull_alignment_crick[start+1] = mismatch
    return [succesfull_alignment_wattson, succesfull_alignment_crick]


def genomic_energy_profile(matrix: list[dict], genome: str, matrix_type: str):
    energy_profile_wattson = []
    energy_profile_crick = []
    genome_to_reverse = genome
    reversed_genome = ''.join(str(x) for x in reversed(genome_to_reverse))
    for pos in range(len(genome)-len(matrix)):
        sequence = genome[pos:pos+len(matrix)]
        reversed_sequence = reversed_genome[pos:pos+len(matrix)]
        if matrix_type == "ppm":
            score = pwm.score_ppm(
                sequence=sequence, PPM=matrix)[1]
            max_score = pwm.max_score_ppm(matrix)
            energy_profile_wattson.append(score/max_score)
            energy_profile_crick.append(pwm.score_ppm(
                sequence=reversed_sequence, PPM=matrix)[1])
        elif matrix_type == "pwm":
            score = pwm.score_pwm(sequence=sequence, PWM=matrix)[1]
            max_score = pwm.max_score_pwm(matrix)
            energy_profile_wattson.append(score/max_score)
            energy_profile_crick.append(pwm.score_pwm(
                sequence=reversed_sequence, PWM=matrix)[1])
        else:  # if the matrix type is not specified
            print("Matrix type not specified")
            return None

    fig, ax = plt.subplots()
    fig, bx = plt.subplots()
    ax.plot(energy_profile_wattson, label="Energy profile Wattson")
    bx.plot(energy_profile_crick, label="Energy profile Crick")

    # Add title and labels
    ax.set_title('Energy Profile Wattson')
    ax.set_xlabel('Position')
    ax.set_ylabel('Energy')

    bx.set_title('Energy Profile Crick')
    bx.set_xlabel('Position')
    bx.set_ylabel('Energy')

    plt.show()

    return energy_profile_wattson, energy_profile_crick


if __name__ == "__main__":
    print("This script is not meant to be run directly")
    sys.exit(1)
