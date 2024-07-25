
def get_benchmark(path):
    with open(path, "r") as file:
        file_list = [line.strip() for line in file.readlines()]
        benchmark = []
        for line in file_list[6:36]:
            new_line = line.split(sep="\t")
            if new_line[-3] == "+":
                chain = new_line[-6]
            else:
                reversed_nts = ""
                for nt in new_line[-6][::-1]:
                    if nt == "A":
                        reversed_nts += ("T")
                    elif nt == "T":
                        reversed_nts += ("A")
                    elif nt == "C":
                        reversed_nts += ("G")
                    elif nt == "G":
                        reversed_nts += ("C")
                chain = reversed_nts

            if new_line[-1] == "-":
                position = int(new_line[6].replace(',', '')) - int(new_line[3].replace(',', ''))
            else:
                position = int(new_line[3].replace(',', '')) - int(new_line[6].replace(',', ''))

            if position < 0:
                position += 299
                new_chain = ""
                for i in range(len(chain)):
                    new_chain += chain[(i + 1) * -1]
            filtered_line = [new_line[1], new_line[2],
                             position,
                             chain]
            benchmark.append(filtered_line)
    return benchmark


if __name__ == "__main__":
    benchmark = get_benchmark("_input_files/benchmark_marboxes.txt")
    print(benchmark)
