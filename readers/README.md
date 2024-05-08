# Readers Module

This module contains functions to read files in different formats, so that they can be used as input for other functions in the package.

# `fasta_reader.py`


## `fasta_reader` Function

This function reads a fasta file and returns each sequence as a string in a list.

### Parameters

- `fasta_file_path` (str): Path to the fasta file.
- `raw` (bool): If True, reads every line as a new sequence. Used for .seq files. Default is False.

### Returns

- `list[str]`: List of nucleotide sequences in the fasta file.

### Usage

The `fasta_reader` function can be used to read a fasta file and return each sequence as a string in a list. The function takes as input the path to the fasta file.

Here's an example of how to use the `fasta_reader` function:

```python
fasta_file_path = "path_to_your_fasta_file.fasta"
sequences = fasta_reader(fasta_file_path)
```

In this example, `sequences` will be a list of strings, where each string is a nucleotide sequence from the fasta file.

# `matrix_reader.py`

This file contains the following function:

## `read_matrix` in `matrix_reader.py`

This function opens and reads a matrix file, saving it as a pandas DataFrame. It ignores lines that start with 'PO' or '#'. The first row of the matrix is used as the column headers in the DataFrame.

### Parameters

- `matrix_file_path` (str): Path to the matrix file.

### Returns

- `pd.DataFrame`: DataFrame with the matrix.

### Usage

Here's an example of how to use the `read_matrix` function:

```python
import matrix_reader

df = matrix_reader.read_matrix("path_to_your_matrix_file")
```

This will return a pandas DataFrame with the matrix from the file at "path_to_your_matrix_file".