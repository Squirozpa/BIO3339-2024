""" 
Script designed to recieve any length nucleotidic data, and create a file of PWM
"""

import sys
import os


def open_file(file_path):
    """Function that recieves a path to open file and return it as a list"""

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
    """Function recieves the list of sequences and creates a list of dictionaries with the length of the first string, then returns the count of nts for each position"""

    nts_distribution = [{'A': 0, 'C': 0, 'G': 0, 'T': 0}
                        for n in range(len(nts_list[0]))]
    for line in nts_list:
        # for each "sequence" (line)
        for i in range(len(line)):
            # for every position, adds a count to the corresponding nts of that position
            nts_distribution[i][line[i]] += 1
    return nts_distribution


def pwmer(nts_list):
    """Function that recieves the list already counted and asks for vertical or horizontal PWM, then returns the PWM in the specified format"""

    marker = input("Please select vertical or horizontal (H/V) : ")

    if marker.lower() == "h":
        nts_headers = ["A", "C", "G", "T"]
        final = f"PO    {'  '.join(map(str, range(1, len(nts_list) + 1)))}\n"
        for nts in range(4):
            # to traverse 4 times one for each nts
            nts_check = nts_headers.pop(0)
            line = f"{nts_check}    "
            for position in range(len(nts_list)):
                line += f"{str(nts_list[position][nts_check])}    "
            final += f"{line}\n"

    elif marker.lower() == "v":
        final = "PO    A    C    G   T\n"
        for posicion in range(len(nts_list)):
            line = f"{str(posicion + 1)}    {nts_list[posicion]['A']}    {nts_list[posicion]['C']}    {nts_list[posicion]['G']}    {nts_list[posicion]['T']}"
            final += f"{line}\n"

    else:
        print("Please use a single letter V/H for vertical or horizontal")
        return pwmer(nts_list)
    return final


def saver(pwm_str):
    """Function that recieves tha PWM str and saves it as a file with the specified name"""

    file_name = input("Please type the file name to be save as: \n")
    file = open(f"{file_name}.txt", "w")
    file.write(pwm_str)
    return file_name


if __name__ == "__main__":
    path_folder = os.path.dirname(__file__)
    file_name = sys.argv[1]
    path = os.path.join(path_folder, file_name)

    file = open_file(path)
    counted = counter(file)
    pwmed = pwmer(counted)
    name = saver(pwmed)
    print(f"File {name}.txt succesfully created")
