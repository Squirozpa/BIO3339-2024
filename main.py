"""

"""

# Standard Library Imports
import pandas as pd
import argparse
# Local Library Imports
import writers.energy_writer as ew
import writers.alignment_writer as aw
import writers.matrix_writer as mw
import writers.fasta_writer as fw
import matrix.matrix as mm
import matrix.scores as ms
import readers.fasta_reader as fr
import readers.matrix_reader as mr
import alignment.alignment as al
import consensus.consensus as cs
import energy.energy as en
import energy.graphs as eg
####################################################################################################

pd.set_option('future.no_silent_downcasting', True)
pd.set_option('display.float_format', '{:.2f}'.format)


def prepend_input_path(file_path):
    if file_path.startswith('='):
        input_file_path = f'_output_files/{file_path[1:]}'
    else:
        input_file_path = f'_input_files/{file_path}'
    return input_file_path


def main():
    parser = argparse.ArgumentParser(description='Process some files.')
    parser.add_argument('command', type=str, choices=[
                        'matrix', 'consensus', 'alignment', 'energy'], help='The command to run')
    parser.add_argument('--input_file', type=str, required=True,
                        help='Path to the fasta file with the input_file')
    parser.add_argument('--output_name', type=str,
                        required=True, help='Base name for the output files')
    parser.add_argument('--vertical', action=argparse.BooleanOptionalAction, default=False,
                        help='If the matrix should be vertical or horizontal')
    parser.add_argument('--gc_content', action=argparse.BooleanOptionalAction, default=False,
                        help='If the matrix should be normalized by the GC content')
    parser.add_argument('--raw', action=argparse.BooleanOptionalAction, default=False,
                        help='If the MSA is a .seq file')
    parser.add_argument('--threshold', type=float,
                        default='0.5', help='threshold for the consensus')
    parser.add_argument('--prioritize', action=argparse.BooleanOptionalAction, default=False,
                        help='If the consensus should prioritize upper thresholds')
    parser.add_argument('--matrix_type', type=str, choices=[
                        "weight", "relative"], default='weight', help='Type of matrix to use for the consensus')
    parser.add_argument('--target_file', type=str,
                        help='Path to the fasta file with the target sequences for the alignment')
    parser.add_argument('--max_mismatch', type=int, default=None,
                        help='Maximum number of mismatches for the alignment')
    parser.add_argument('--normalized', type=str, choices=[
                        "normalized", "relative", "none"], default='normalized', help='If the energy values should be normalized or relative')
    parser.add_argument('--cut_off', type=float, default=None,
                        help='Cut off value for the energy')
    parser.add_argument('--reversed', action=argparse.BooleanOptionalAction, default=True,
                        help='If the energy values should be calculated for the reverse sequence')
    parser.add_argument('--only_reversed', action=argparse.BooleanOptionalAction, default=False,
                        help='If only the reverse energy values should be calculated')
    parser.add_argument('--graph', action=argparse.BooleanOptionalAction, default=True,
                        help='If the energy graph should be generated')
    parser.add_argument('--graph_length', type=int, default='400',
                        help='Length of the graph for the energy profile')

    args = parser.parse_args()
    input_file_path = prepend_input_path(args.input_file)
    # Prepend standard paths
    output_name = f'_output_files/{args.output_name}'

    if args.command == 'matrix':
        input_file = fr.read_fasta(input_file_path, raw=args.raw)
        # Generate matrices
        absolute_matrix = mm.absolute_matrix(input_file)
        relative_matrix = mm.relative_matrix(input_file)
        if args.gc_content:
            target_file_path = prepend_input_path(args.target_file)
            genome_dict = mm.genomic_count(input_file)
            weight_matrix = mm.weight_matrix(
                input_file, genome_frequencey=genome_dict)
        else:
            weight_matrix = mm.weight_matrix(input_file)

        # Save matrices
        mw.write_matrix(absolute_matrix,
                        f'{output_name}_pfm', vertical=args.vertical)
        mw.write_matrix(relative_matrix, f'{output_name}_ppm',
                        max_value=ms.max_score_relative_matrix(relative_matrix), vertical=args.vertical, decimals=True)
        mw.write_matrix(weight_matrix, f'{output_name}_pwm',
                        max_value=ms.max_score_weight_matrix(weight_matrix), vertical=args.vertical, decimals=True)

    if args.command == 'consensus':
        # Read input_file
        input_file = mr.read_matrix(input_file_path)

        # Generate consensus
        consensus = cs.consensus(input_file, threshold=args.threshold,
                                 prioritize_upper=args.prioritize, matrix_type=args.matrix_type)

        fw.write_fasta(
            consensus, f"consensus_{args.output_name}_threshold:{args.threshold}_matrixtype:{args.matrix_type}", f'{output_name}_consensus')

    if args.command == 'alignment':
        target_file_path = f'_input_files/{args.target_file}'
        target_sequence = fr.read_fasta(target_file_path)
        input_file = fr.read_fasta(input_file_path)
        alignment = al.align_both_directions(
            input_file[0], target_sequence[0], max_mismatch=args.max_mismatch)
        aw.write_alignment(alignment,
                           f'{output_name}_alignment', max_mismatch=args.max_mismatch)

    if args.command == 'energy':
        target_file_path = f'_input_files/{args.target_file}'
        matrix = mr.read_matrix(input_file_path)
        target = fr.read_fasta(target_file_path)
        energy = en.genomic_energy(matrix=matrix, sequence=target[0], matrix_type=args.matrix_type,
                                   threshold=args.cut_off, normalized=args.normalized, reversed=args.reversed, only_reversed=args.only_reversed)
        ew.write_energy(energy, f'{output_name}_energy.tsv')
        if args.graph:
            eg.plot_energy_profile(energy, output_name,
                                   graph_length=args.graph_length, length_matrix=len(matrix))


if __name__ == "__main__":
    main()
