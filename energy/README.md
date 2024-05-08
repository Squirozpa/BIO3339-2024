# Energy Module

This module provides functions for calculating the energy of a sequence in a genomic sequence, given a matrix, and plotting the energy profile of a sequence.

# `energy.py`

This file contains the following functions:

## `genomic_energy` in `energy.py`

This function calculates the energy of a sequence in a genomic sequence, given a matrix.

### Parameters

- `matrix` (pd.DataFrame): The weight or relative matrix used for energy calculation.
- `sequence` (str): The genomic sequence for which energy needs to be calculated.
- `matrix_type` (str, optional): The type of matrix used for energy calculation. Default is "weight".
- `threshold` (int, optional): Threshold for appending values bigger than the threshold. Default is None.
- `normalized` (str, optional): Flag indicating whether the energy values should be normalized. Options are "normalized", "relative", or "none". Default is "normalized".
- `reversed` (bool, optional): Flag indicating whether the reversed energy values should also be calculated. Default is True.

### Returns

- pd.DataFrame: A DataFrame of energy values for the given sequence. If `reversed` is True, it also returns reversed energy values.

### Raises

- ValueError: If the `matrix_type` is not "weight" or "relative".

### Usage

Here's an example of how to use the `genomic_energy` function:

```python
matrix = ...  # Assume this is a pandas DataFrame with the weight or relative matrix
sequence = "ATCG"
genomic_energy = energy.genomic_energy(matrix, sequence)
```

# `graphs.py`

This file contains the following functions:

## `plot_energy_profile` in `graphs.py`

This function plots the energy profile of a sequence.

### Parameters

- `energy_profile` (pd.DataFrame): A pandas DataFrame with the energy profile.
- `filename` (str): The filename to save the plot.
- `title` (str, optional): The title of the plot. Default is "Energy Profile".
- `graph_length` (int, optional): The length of the graph. Default is 400.

### Usage

Here's an example of how to use the `plot_energy_profile` function:

```python
energy_profile = ...  # Assume this is a pandas DataFrame with the energy profile
filename = "energy_profile.png"
title = "Energy Profile of Sequence"
graph_length = 500
graphs.plot_energy_profile(energy_profile, filename, title, graph_length)
```

This function will plot the energy profile of the sequence and save the plot to a file. The plot will also be displayed.
