# Writers Module

This module provides functions for writing various data structures to files. The writers module includes functions for writing sequences, matrices, alignments, and energy profiles to files in different formats.

# `energy_writer.py` 

## `write_alignment` Function

This function writes the alignment of sequences to a file and saves it as a `.txt` file.

### Parameters

- `alignment` (list[dict]): The alignment as a list of dictionaries. Each dictionary represents an alignment and contains the mismatches for each alignment.
- `output_file` (str): The path to the output file.
- `max_mismatch` (int, optional): The maximum number of mismatches allowed. If provided, this value is added to the output file.

### Returns

- `output_file` (str): The path to the output file.

### Usage

The `write_alignment` function can be used to write the alignment of sequences to a file. The function takes as input an alignment (a list of dictionaries), the path to the output file, and optionally the maximum number of mismatches allowed.

Here's an example of how to use the `write_alignment` function:

```python
alignment = [[{1: [0, 2], 0: [4]}, {2: [0, 1], 0: [3]}]]
output_file = "alignment.txt"
max_mismatch = 2

write_alignment(alignment, output_file, max_mismatch)
```

### Output_file

```txt
Alignment Alignment Wattson::
Mismatches(1): [0,2]
Mismatches(0): [4]

Alignment Alignment Crick::
Mismatches(2): [0,1]
Mismatches(0): [3]

Max Mismatches: 2
```

# energy_writer.py

## `write_energy` Function

This function exports the genomic energy profile to a file and saves it as a `.tsv` file.

### Parameters

- `genomic_energy` (pd.DataFrame): The DataFrame of energy values for the genomic sequence.
- `output_file` (str): The path to the output file.

### Returns

- `output_file` (str): The path to the output file.

### Usage

The `write_energy` function can be used to export the genomic energy profile to a file. The function takes as input a DataFrame of energy values for the genomic sequence and the path to the output file.

Here's an example of how to use the `write_energy` function:

```python
genomic_energy = genomic_energy(matrix, sequence, matrix_type)
output_file = "energy"

write_energy(genomic_energy, output_file)
```

### Output_file

```tsv
energy  reversed_energy
0.3   0.1
0.1   0.2
0.4   0.1
0.9   0
1.1   2
```

# `fasta_writer.py`

## `write_fasta` Function

This function exports a sequence in FASTA format and saves it as a `.fasta` file.

### Parameters

- `sequence` (str): Sequence to be exported.
- `sequence_id` (str): ID of the sequence.
- `output_file` (str): Name of the output file.

### Returns

- `output_file` (str): Name of the output file.

### Usage

The `write_fasta` function can be used to export a sequence in FASTA format. The function takes as input a sequence, the ID of the sequence, and the name of the output file.

Here's an example of how to use the `write_fasta` function:

```python
sequence = "ATCGTTGACTGATCGTACGATCGTACG"
sequence_id = "seq1"
output_file = "sequence"

write_fasta(sequence, sequence_id, output_file)
```

### Output_file

```txt
>seq1
ATCGTTGACTGATCGTACGATCGTACG
```

# `matrix_writer.py`

## `write_matrix` Function

This function exports a matrix in tab-separated format and saves it as a `.tsv` file.

### Parameters

- `matrix` (pd.DataFrame): Matrix to be exported.
- `output_file` (str): Path to save the output file.
- `max_value` (float, optional): Maximum value to be added at the end of the file. Defaults to None.
- `vertical` (bool, optional): Flag to indicate if the matrix should be exported vertically. Defaults to False.

### Returns

- `output_file` (str): Name of the output file.

### Usage

The `write_matrix` function can be used to export a matrix in tab-separated format. The function takes as input a matrix, the path to save the output file, an optional maximum value to be added at the end of the file, and an optional flag to indicate if the matrix should be exported vertically.

Here's an example of how to use the `write_matrix` function:

```python
import pandas as pd

matrix = pd.DataFrame({
    'A': [1, 2, 3],
    'C': [4, 5, 6],
    'G': [7, 8, 9],
    'T': [3, 1, 5]
})

output_file = "matrix"
max_value = 9
vertical = True

write_matrix(matrix, output_file, max_value, vertical)
```

### Output_file

```txt
PO	A	C	G	T
0	1	4	7	3
1	2	5	8	1
2	3	6	9	5
#max value: 9
```
