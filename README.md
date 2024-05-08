# Transcription Factor Binding Site Analysis

This project is a command-line tool for processing fasta files and performing various operations such as matrix generation, consensus, alignment, and energy calculation.

## Getting Started

Download the project and its prerequisites. Follow the instructions in usage to run the program.

### Prerequisites

You need to have the following installed on your machine:

`pandas`

``` bash
pip install pandas
```

`matplotlib`

``` bash
pip install matplotlib
```

`numpy`

``` bash
pip install numpy
```

Note: matplotlib auto installs numpy as a dependency.

## Usage

Here is how to use the `main.py` script:

```bash
python main.py <command> --input_file <input_file> --output_name <output_name> [OPTIONS...]
```

Where:

```txt
<command> is one of 'matrix', 'consensus', 'alignment', 'energy'.
<input_file> is the path to the fasta file with the input sequences.
<output_name> is the base name for the output files.
[OPTIONS...] are additional options, such as:
--vertical: If the matrix should be vertical (default is False).
--gc_content: If the matrix should be normalized by the GC content (default is False). Note: This option requeries a target file.
--threshold: Threshold for the consensus (default is 0.5).
--prioritize: If the consensus should prioritize upper thresholds (default is False).
--matrix_type: Type of matrix to use for the consensus or energy calculations (default is 'weight'). Options are 'weight', 'relative'.
--target_file: Path to the fasta file with the target sequences for the alignment or energy calculations.
--max_mismatch: Maximum number of mismatches for the alignment (default is None).
--normalized: If the energy values should be normalized or relative (default is 'normalized'). Options are 'normalized', 'relative', 'none'.
--cut_off: Cut off value for the energy (default is None).
--no-reversed: If the energy values should not be calculated for the reverse sequence (default is True).
--no-graph: If the energy graph should not be generated (default is True).
--graph_length: Length of the graph for the energy profile (default is 400).
```

Input files should be in the folder `_input_files` and output files will be saved in the folder `_output_files`.  
All output files have automatically their extension added according to the command used. Dont add the extension to the `--output_name`. Else it will save as `output_name.extension.extension`.  

Note: Using "=" before an `--input_file`, selects that input file from the `_output_files` folder. Makes it easier to use the next commands without moving the output to the input folder. Example: `--input_file output_name.fasta`.

## Functions

### `Matrix`

The command `matrix` generates a matrix from the input sequences. The matrix can be horizontal or vertical and normalized by the GC content.  
Note: For all the next calculations matrix are used as vertical. Matrix inputs can be horizontal or vertical, but if the matrix is horizontal, it will be transposed. If the horizontal matrix has less length than 5, it wont be transposed, and the input matrix must be manually set to vertical.

#### Arguments for `matrix`

`--vertical`
If the matrix should be vertical or horizontal (default is False).  
Options: are True or False.

`--gc_content`
If the weight matrix should be normalized by the GC content (default is False).
Options are True or False.  
Note: This option requeries a target file. With the option `--target_file` the target sequences are used to calculate the GC content.

`--raw`
If the input file is a raw file (default is False).
Options are True or False.
Used if the input file is a .seq file.

### `Consensus`

The command `consensus` generates a consensus sequence from the input sequences following IUPAC codes for degeneracy. The consensus can be generated with a threshold and prioritize upper thresholds. The consensus can be generated using a weight matrix or a relative matrix. The type of matrix used changes the way the consensus is calculated.  

#### The calculations for relative matrix are

- It calculates 3 different "__adjusted thresholds__". Each representing a different amount of nucleotides to be considered per position. Using the following formulae, where k is the threshold set as the input for the command:

  - upper_threshold: 0.25 + (0.75 * k)
  - middle_threshold: 0.25 + (0.25 * k)
  - lower_threshold: 0.25 + (0.08 * k)

  Each __adjusted threshold__ represents the amount of nucleotides that must be present in the position, times the ratio of the threshold.
  For k = 1, the __adjusted thresholds__ are 1, 0.5, 0.33.
  Meaning only the maximum amount of nucleotides must be present in the position to be considered.  
  Note: More than one statement may be true for the same position when using k =< 0.33. So thats why there is a `--prioritize_upper` option.

- Then, using the relative matrix, it selects the nucleotides that are present in the position according to the threshold set. If theres only one that exceeds __upper_threshold__, it is selected. If there are 2 that exceed __middle_threshold__, they are selected. If there are 3 that exceed __lower_threshold__, they are selected. The order of selection is reversed if the `--prioritize` option is set to False.

- Selected nucleotides are converted to the IUPAC code. If there are no nucleotides that exceed the thresholds, the consensus is set to 'N'.

- The consensus is saved as fasta format.

Note: Using the relative matrix is not very usefull because of equiprobability bias, but it is a good way to see the difference between the weight matrix and the relative matrix.

#### The calulations for weight matrix are

- It calculates 3 different "__adjusted thresholds__". Each representing a different amount of nucleotides to be considered per position. Using the following formulae, where k is the threshold set as the input for the command:
  - upper_threshold: k*log2(1/0.25) = k*(log2(4)) = 2k
  - middle_threshold: k*log2(0.5/0.25) = k*(log2(2)) = k
  - lower_threshold: k*log2(0.33/0.25) = k*(log2(1.33)) ≈ 0.4k

Each __adjusted threshold__ represents the "information content" per position for each nucleotide to be selected.
For k = 1, the __adjusted thresholds__ are 2, 1, 0.4.
Meaning only the maximum amount of "information content" must be present in the position to be considered.  
Note: More than one statement may be true for the same position when using k =< 0.5. So thats why there is a `--prioritize_upper` option.

- Then, using the weight matrix, it selects the nucleotides that are present in the position according to the threshold set. If theres only one that exceeds __upper_threshold__, it is selected. If there are 2 that exceed __middle_threshold__, they are selected. If there are 3 that exceed __lower_threshold__, they are selected. The order of selection is reversed if the `--prioritize` option is set to False.

- Selected nucleotides are converted to the IUPAC code. If there are no nucleotides that exceed the thresholds, the consensus is set to 'N'.

- The consensus is saved as fasta format.

Note: Using the weight matrix is better, since it uses the information content of the matrix to select the nucleotides. Meaning it accounts for the raririty of the nucleotides in the alignment. The thresholds are now more meaningful, since they are scaled based on the information content of the matrix. The problem of multiple statements being true for the same and its even more probable to happen, but since the information content is scaled, the probability of this happening is lower, and the thresholds are now more sensible at higher values of k, meaning that it is less necessary to lower the k value to see results.

#### Arguments for `consensus`

`--threshold`
Threshold for the consensus (default is 0.5). This is the k value specified above.

`--prioritize`
If the consensus should prioritize upper thresholds (default is False).
Options are True or False.  
Note: Only necessary when using k =< 0.5 for weight matrix or k =<0.33 when using relative. Setting it to True means that the upper thresholds are prioritized, meaning that the lower thresholds are only considered if the upper thresholds are not met. So results nucleotides are less ambigous.

`--matrix_type`
Type of matrix to use for the consensus or energy calculations (default is 'weight').
Options are 'weight', 'relative'.

### `Alignment`

The command `alignment` aligns the input sequences with the target sequences. The alignment is performed using simple alignment algorithm. The alignment can have a maximum number of mismatches. Algorithm doesnt account for gaps. Can interpet ambigous codes, for query sequence The alignment is saved as a txt file.

#### Arguments for `alignment`

`--target_file`
Path to the fasta file with the target sequences for the alignment.

`--max_mismatch`
Maximum number of mismatches for the alignment (default is None).

### `Energy`

The command `energy` calculates the energy of the input matrix with the target sequences. The energy can be normalized or relative more info below. The energy can have a cut off value meaning values scores of the matrix are not selected (they are not the normalized or relativized values, option coming soon). The energy can be calculated for the reverse sequence. The energy is saved as a tsv file. And the option to plot to matplotlib is available.

#### The calculations for energy are

This function calculates the energy per position of the target sequence using the input matrix. The energy is calculated as the sum (for weight matrix) or product (for relative matrix) of the values of the matrix for the nucleotides in the position.  
__Normalization__ means the score per position is divided by the maximum score of the matrix.  
__Relative__ scores are calculated as follows:

$$
score_r = \frac{score - min}{max - min}
$$

Where min and max are the minimum and maximum values of the matrix.

Cut off value is applied before normalization or relativization.
Meaning that the values that are below the cut off value are not considered for the normalization or relativization. Default is None, meaning that all values are considered. Recommended option is to use 0 so that values that are below 0 are not considered in a weight matrix.

#### Arguments for `energy`

`--normalized`
If the energy values should be normalized or relative (default is 'normalized').
Options are 'normalized', 'relative', 'none'.

`--cut_off`
Cut off value for the energy (default is None).

`--reversed`
Flag to calculate the energy for the reverse sequence (default is True).

`--graph`
If the energy graph should be generated (default is True).
Options are True or False.  
Takes more time to plot the graph, so if you dont want the graph, set it to False.

`--graph_length`
Length of the graph for the energy profile (default is 400).

## Disclaimer

This project is a university assignment and is not intended for distribution or use outside of the course. The code is provided as-is and may contain bugs or errors. Use at your own risk.

## Authors

__Sebastián Quiroz__  
Undergraduate student at Pontificia Universidad Católica de Chile  
Email: <squirozpa@uc.cl>

__Francisco Melo__  
Dr. in Structural Biology  
Professor at Pontificia Universidad Católica de Chile  
Email: <fmelo@bio.puc.cl>
