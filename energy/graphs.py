"""
Functions to plot the energy profile of a sequence
"""

# Standard Library Imports
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
# Local Library Imports
import energy.energy as en
####################################################################################################


def plot_energy_profile(energy_profile: pd.DataFrame, filename: str, title: str = "Energy Profile", graph_length: int = 400, length_matrix: int = 0):
    """Function that plots the energy profile of a sequence

    Args:
        energy_profile (pd.DataFrame): A pandas DataFrame with the energy profile
        title (str): The title of the plot
        filename (str): The filename to save the plot
    """
    # Check if there's an 'energy' column in the DataFrame
    if 'energy' in energy_profile.columns:
        # Add null values to the end of the energy column
        energy = energy_profile['energy']._append(
            pd.Series([np.nan]*length_matrix)).reset_index(drop=True)
        energy.name = 'energy'
    else:
        energy = None

    # Check if there's a 'reversed_energy' column in the DataFrame
    if 'reversed_energy' in energy_profile.columns:
        # Add null values to the start of the reversed_energy column
        reversed_energy = pd.Series(
            [np.nan]*length_matrix)._append(energy_profile['reversed_energy']).reset_index(drop=True)
        reversed_energy.name = 'reversed_energy'
    else:
        reversed_energy = None

    # Concatenate the series along axis=1 and fill the null values
    if energy is not None and reversed_energy is not None:
        energy_profile = pd.concat([energy, reversed_energy], axis=1)
    elif energy is not None:
        energy_profile = energy.to_frame()
    elif reversed_energy is not None:
        energy_profile = reversed_energy.to_frame()

    # Determine the number of subplots
    n = len(energy_profile) // graph_length
    if len(energy_profile) % graph_length != 0:
        n += 1
    # Create a figure and axes
    fig, axs = plt.subplots(n, figsize=(19, 9.2))
    # Determine the global minimum and maximum y values
    if "reversed_energy" in energy_profile.columns and "energy" in energy_profile.columns:
        global_min = -energy_profile['reversed_energy'].max()
        global_max = energy_profile['energy'].max()
    elif "reversed_energy" in energy_profile.columns:  # If there is no energy column
        global_min = energy_profile['reversed_energy'].min()
        global_max = energy_profile['reversed_energy'].max()
    elif "energy" in energy_profile.columns:  # If there is no reversed_energy column
        global_min = energy_profile['energy'].min()
        global_max = energy_profile['energy'].max()
    # Plot the energy profile in each subplot
    for i in range(n):
        start = i * graph_length
        # If this is the last iteration, set end to the length of the DataFrame
        if i == n - 1:
            end = len(energy_profile)
        else:
            end = start + graph_length
        x_values = range(start, end)
        if 'energy' in energy_profile.columns:
            axs[i].plot(x_values, energy_profile['energy'][start:end])
        if 'reversed_energy' in energy_profile.columns and "energy" in energy_profile.columns:
            axs[i].plot(x_values, -energy_profile['reversed_energy']
                        [start:end], color='red')
        elif 'reversed_energy' in energy_profile.columns:
            axs[i].plot(x_values, energy_profile['reversed_energy']
                        [start:end], color='red')

        axs[i].set_title(f'{title} {i*graph_length} - {end}')
        axs[i].set_xlabel('Position')
        axs[i].set_ylabel('Energy')
        axs[i].set_ylim(global_min, global_max)
        # Set the x-axis limit to graph_length for all graphs
        axs[i].set_xlim(start, start + graph_length)

    # Save the plot to a file
    plt.tight_layout()
    plt.savefig(filename)
    # Show the plot
    plt.show()
