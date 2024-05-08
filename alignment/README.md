# Alignment Module

This module provides functions for aligning sequences.

# `alignment.py`

This file contains the following function:

## `align` in `alignment.py`

This function aligns a query sequence with a target sequence and returns a dictionary of alignments.

### Parameters

- `query_sequence` (str): The query sequence to align.
- `target_sequence` (str): The target sequence to align against.
- `max_mismatch` (int, optional): The maximum number of allowed mismatches. Defaults to 0.

### Returns

- `dict[int, list]`: A dictionary where the keys represent the number of mismatches and the values are lists of start positions.

### Usage

Here's an example of how to use the `align` function:

```python
query_sequence = "ATCG"
target_sequence = "GCTA"
max_mismatch = 1
alignments = align(query_sequence, target_sequence, max_mismatch)
```

This function will align the query sequence with the target sequence and return a dictionary of alignments. The keys of the dictionary represent the number of mismatches and the values are lists of start positions.  
Notes:

- The start positions are 0-based.
- For reversed alignment, the output will be in the starting position of the 5'-3' direction of the target sequence. So if the target sequence is "AACC" and the query sequence is "GGTT", the output will be 3 (starting position).

```txt
{0: [3]}
```
