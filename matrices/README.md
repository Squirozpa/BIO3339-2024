# Matrix

Functions to generate and manipulate matrices.

# `matrix.py`

This file contains the functions to generate matrices from MSA.

## `genomic_count` in `matrix.py`

This function calculates the frequency of each nucleotide in a DNA sequence and its reverse complement.

### Parameters

- `sequence` (str): DNA sequence.

### Returns

- `dict`: Dictionary with the frequency of each nucleotide.

### Usage

Here's an example of how to use the `genomic_count` function:

```python
frequency = matrix.genomic_count("ATCG")
```

This will return a dictionary with the frequency of each nucleotide in the sequence "ATCG" and its reverse complement.

## `absolute_matrix` in `matrix.py`

This function generates an absolute matrix from a multiple sequence alignment.

### Parameters

- `alignment` (list): A list of sequences, previously aligned.

### Returns

- `pd.DataFrame`: DataFrame with the absolute matrix.

### Usage

Here's an example of how to use the `absolute_matrix` function:

```python
alignment = ["ATCG", "GCTA", "CGAT", "TACG"]
abs_matrix = matrix.absolute_matrix(alignment)
```

This will return a pandas DataFrame with the absolute matrix from the alignment.

## `relative_matrix` in `matrix.py`

This function generates a relative matrix from a multiple sequence alignment.

### Parameters

- `alignment` (list): A list of sequences, previously aligned.

### Returns

- `pd.DataFrame`: DataFrame with the relative matrix.

### Usage

Here's an example of how to use the `relative_matrix` function:

```python
alignment = ["ATCG", "GCTA", "CGAT", "TACG"]
rel_matrix = matrix.relative_matrix(alignment)
```

This will return a pandas DataFrame with the relative matrix from the alignment.

## `entropy_matrix` in `matrix.py`

This function generates an entropy matrix from a multiple sequence alignment.

### Parameters

- `alignment` (list): A list of sequences, previously aligned.

### Returns

- `pd.DataFrame`: DataFrame with the entropy matrix.

### Usage

Here's an example of how to use the `entropy_matrix` function:

```python
alignment = ["ATCG", "GCTA", "CGAT", "TACG"]
entropy_matrix = matrix.entropy_matrix(alignment)
```

This will return a pandas DataFrame with the entropy matrix from the alignment.

## `weight_matrix` in `matrix.py`

This function generates a weight matrix from a multiple sequence alignment.

### Parameters

- `alignment` (list): A list of sequences, previously aligned.
- `genome_frequency` (dict, optional): A dictionary with the frequency of each nucleotide. Defaults to {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}.

### Returns

- `pd.DataFrame`: DataFrame with the weight matrix.

### Usage

Here's an example of how to use the `weight_matrix` function:

```python
alignment = ["ATCG", "GCTA", "CGAT", "TACG"]
weight_matrix = matrix.weight_matrix(alignment)
```

This will return a pandas DataFrame with the weight matrix from the alignment.

# `scores.py`

This file contains the following functions:

## `score_relative_matrix` in `scores.py`

This function calculates the score of a sequence given a relative matrix.

### Parameters

- `sequence` (str): DNA sequence.
- `matrix` (pd.DataFrame): A pandas DataFrame with the relative matrix.

### Returns

- `float`: The score of the sequence given a matrix.

### Usage

Here's an example of how to use the `score_relative_matrix` function:

```python
sequence = "ATCG"
rel_matrix = ...  # Assume this is a pandas DataFrame with the relative matrix
score = scores.score_relative_matrix(sequence, rel_matrix)
```

## `score_weight_matrix` in `scores.py`

This function calculates the score of a sequence given a weight matrix.

### Parameters

- `sequence` (str): DNA sequence.
- `matrix` (pd.DataFrame): A pandas DataFrame with the weight matrix.

### Returns

- `float`: The score of the sequence given a matrix.

### Usage

Here's an example of how to use the `score_weight_matrix` function:

```python
sequence = "ATCG"
weight_matrix = ...  # Assume this is a pandas DataFrame with the weight matrix
score = scores.score_weight_matrix(sequence, weight_matrix)
```

## `max_score_relative_matrix` in `scores.py`

This function calculates the maximum score of a relative matrix.

### Parameters

- `matrix` (pd.DataFrame): A pandas DataFrame with the relative matrix.

### Returns

- `float`: The maximum score of the matrix.

### Usage

Here's an example of how to use the `max_score_relative_matrix` function:

```python
import scores

rel_matrix = ...  # Assume this is a pandas DataFrame with the relative matrix
max_score = scores.max_score_relative_matrix(rel_matrix)
```

## `max_score_weight_matrix` in `scores.py`

This function calculates the maximum score of a weight matrix.

### Parameters

- `matrix` (pd.DataFrame): A pandas DataFrame with the weight matrix.

### Returns

- `float`: The maximum score of the matrix.

### Usage

Here's an example of how to use the `max_score_weight_matrix` function:

```python
weight_matrix = ...  # Assume this is a pandas DataFrame with the weight matrix
max_score = scores.max_score_weight_matrix(weight_matrix)
```

## `min_score_relative_matrix` in `scores.py`

This function calculates the minimum score of a relative matrix.

### Parameters

- `matrix` (pd.DataFrame): A pandas DataFrame with the relative matrix.

### Returns

- `float`: The minimum score of the matrix.

### Usage

Here's an example of how to use the `min_score_relative_matrix` function:

```python
rel_matrix = ...  # Assume this is a pandas DataFrame with the relative matrix
min_score = scores.min_score_relative_matrix(rel_matrix)
```

## `min_score_weight_matrix` in `scores.py`

This function calculates the minimum score of a weight matrix.

### Parameters

- `matrix` (pd.DataFrame): A pandas DataFrame with the weight matrix.

### Returns

- `float`: The minimum score of the matrix.

### Usage

Here's an example of how to use the `min_score_weight_matrix` function:

```python
weight_matrix = ...  # Assume this is a pandas DataFrame with the weight matrix
min_score = scores.min_score_weight_matrix(weight_matrix)
```

## `relative_score` in `scores.py`

This function calculates the relative score of a sequence given the minimum and maximum scores.

### Parameters

- `min_score` (float): The minimum score of the matrix.
- `max_score` (float): The maximum score of the matrix.
- `score` (float): The score of the sequence.

### Returns

- `float`: The relative score of the sequence.

### Usage

Here's an example of how to use the `relative_score` function:

```python
min_score = ...  # Assume this is the minimum score of the matrix
max_score = ...  # Assume this is the maximum score of the matrix
score = ...  # Assume this is the score of the sequence
relative_score = scores.relative_score(min_score, max_score, score)
```
