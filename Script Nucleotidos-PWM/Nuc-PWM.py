""" 
Script designed to recieve any length nucleotidic data, and create a file of PWM
"""

import sys

"Function to open file and save it as a list"


def open_file(file_path):
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


def counter(nts_list):
    nts_distribution = [{'A': 0, 'T': 0, 'G': 0, 'C': 0}
                        for n in range(len(nts_list[0]))]
    for line in nts_list:
        for i in range(len(line)):
            nts_distribution[i][line[i]] += 1
    return nts_distribution


if __name__ == "__main__":
    path = "Script Nucleotidos-PWM\marboxes_27nts.fasta"
    caca = open_file(path)
    poto = counter(caca)
    print((poto))
