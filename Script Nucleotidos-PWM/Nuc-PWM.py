""" 
Script designed to recieve any length nucleotidic data, and create a file of PWM
"""

import sys

"Function that recieves a path to open file and return it as a list"


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


"Function recieves the list of sequences and creates a list of dictionaries with the length of the first string, then returns the count of nts for each position"


def counter(nts_list):
    nts_distribution = [{'A': 0, 'C': 0, 'G': 0, 'T': 0}
                        for n in range(len(nts_list[0]))]
    for line in nts_list:
        # for each "sequence" (line)
        for i in range(len(line)):
            # for every position, adds a count to the corresponding nts of that position
            nts_distribution[i][line[i]] += 1
    return nts_distribution


"Function that recieves the list already counted and asks for vertical or horizontal PWM, then returns the PWM in the specified format"


def pwmer(nts_list):
    marker = input("Please select vertical or horizontal (H/V) : ")

    if marker.lower() == "v":
        nts_headers = ["A", "C", "G", "T"]
        final = f"      {'  '.join(map(str, range(1, len(nts_list) + 1)))}\n"
        for nts in range(4):
            nts_check = nts_headers.pop(0)
            line = nts_check + "   "
            for position in range(len(nts_list)):
                line += str(nts_list[position][nts_check]) + "   "
            final += f"{line}\n"

    elif marker.lower() == "h":
        final = "   A   C   G   T\n"
        for posicion in range(len(nts_list)):
            line = f"{str(posicion + 1)}    {nts_list[posicion]['A']}    {nts_list[posicion]['C']}    {nts_list[posicion]['G']}    {nts_list[posicion]['T']}"
            final += f"{line}\n"
    return final


if __name__ == "__main__":
    path = "Script Nucleotidos-PWM\marboxes_27nts.fasta"
    caca = open_file(path)
    poto = counter(caca)
    pipi = pwmer(poto)
    print(pipi)
