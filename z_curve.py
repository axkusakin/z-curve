"""This module calculates and plots z-curve (https://doi.org/10.1006/jtbi.1997.0401)"""

import pandas as pd
import plotly.express as px

def delete_description(file):
    """Delete description from FASTA file.
    
    Doesn't work on multifasta.
    """
    total_sequence = []
    with open(file, encoding="utf-8") as sequence:
        lines = sequence.read().split('\n')
    for line in lines:
        if line.find('>') == 0:
            pass
        else:
            total_sequence.append(line)
    return ''.join(total_sequence)

def genome_coordinates(genome):
    """Counts coordinates.
    
    x = R-Y (purine - pyrimidine),
    y = M-K (amino - keto),
    z = W-S (strongH bond - weak-H bond).
    """
    count_x, count_y, count_z = 0, 0, 0
    coordinates = {'x': [], 'y': [], 'z': []}
    for nucleotide in genome:
        if nucleotide == 'A':
            count_x += 1
            count_y += 1
            count_z += 1
        elif nucleotide == 'T':
            count_x -= 1
            count_y -= 1
            count_z += 1
        elif nucleotide == 'G':
            count_x += 1
            count_y -= 1
            count_z -= 1
        elif nucleotide == 'C':
            count_x -= 1
            count_y += 1
            count_z -= 1
        coordinates['x'].append(count_x * 2)
        coordinates['y'].append(count_y * 2)
        coordinates['z'].append(count_z * 2)
    return coordinates

def plot_maker(genome):
    """Plots the z-curve."""
    coordinates = genome_coordinates(genome)
    dataframe_dict = {"x": coordinates['x'],        # Create dataframe from coordinates
                      "y": coordinates['y'],
                      "z": coordinates['z']}
    dataframe = pd.DataFrame.from_dict(dataframe_dict)
    fig = px.line_3d(dataframe, x="x", y="y", z="z")
    return fig.show()


# Дальше код для запуска поделенного генома

def split_genome(genome, sequence_length):
    """Divides the genome into sequences of a given length (sequence_length)
    and averages the parameter (x, y, and z) for each.

    Needed for large genomes, to plot faster.
    ERROR: If genome length is not divisible by sequence_length, the last nucleotides are dropped.
    """
    sequences = []
    for n in range(len(genome) - (len(genome) % sequence_length)):
        if n == 0 or n % sequence_length == 0:
            sequences.append(genome[n : n + sequence_length])
    return sequences

def splitted_genome_coordinates(splitted_genome):
    """Coordinates for the genome divided into sequences.
    
    The function is needed in order to reduce the number of coordinates to construct the z-curve.
    """
    count_x, count_y, count_z = 0, 0, 0
    coordinates = {'x': [], 'y': [], 'z': []}
    for s in range(len(splitted_genome)):
        for n in range(10):
            if splitted_genome[s][n] == 'A':
                count_x += 1
                count_y += 1
                count_z += 1
            elif splitted_genome[s][n] == 'T':
                count_x -= 1
                count_y -= 1
                count_z += 1
            elif splitted_genome[s][n] == 'G':
                count_x += 1
                count_y -= 1
                count_z -= 1
            elif splitted_genome[s][n] == 'C':
                count_x -= 1
                count_y += 1
                count_z -= 1
        coordinates['x'].append(count_x * 2)
        coordinates['y'].append(count_y * 2)
        coordinates['z'].append(count_z * 2)
    return coordinates
