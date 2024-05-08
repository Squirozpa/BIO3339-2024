# Consensus Module

This module provides functions for generating a consensus sequence from a weight or relative matrix.

# `consensus.py`

This file contains the following function:

## `consensus` in `consensus.py`

This function generates a consensus sequence from a weight or relative matrix.

### Parameters

- `matrix` (pd.DataFrame): A pandas DataFrame with the weight or relative matrix.
- `matrix_type` (str, optional): The type of matrix used for consensus sequence generation. Default is "weight".
- `threshhold` (int, optional): The threshold to consider a nucleotide in the consensus sequence. Default is 0.5.
- `prioritize_upper` (bool, optional): Flag indicating whether to prioritize nucleotides that exceed the upper thresholds. Default is False.

### Returns

- `str`: The consensus sequence.

### Usage

Here's an example of how to use the `consensus` function:

```python
matrix = ...  # Assume this is a pandas DataFrame with the weight or relative matrix
threshhold = 0.6
prioritize_upper = True
consensus_sequence = consensus(matrix, "weight", threshhold, prioritize_upper)
```

This function will generate a consensus sequence from the given matrix. The consensus sequence is a string of nucleotides that represents the most common nucleotides at each position in the matrix.

### Output

```python
'NNRAGTCNN'
```
