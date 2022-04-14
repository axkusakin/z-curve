"""This module calculates and plots z-curve (https://doi.org/10.1006/jtbi.1997.0401)"""

import argparse
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
    dataframe_dict = {"x": coordinates['x'],
                      "y": coordinates['y'],
                      "z": coordinates['z']}
    dataframe = pd.DataFrame.from_dict(dataframe_dict)
    fig = px.line_3d(dataframe, x="x", y="y", z="z")
    return fig.show()


# Futher starts code for genome divided into sequences.

def split_genome(genome, sequence_length=10):
    """Divides the genome into sequences of a given length (sequence_length)
    and averages the parameter (x, y, and z) for each.
    Needed for large genomes, to plot faster.
    ERROR: If genome length is not divisible by sequence_length, the last nucleotides are dropped.
    """
    sequences = []
    for nucleotide in range(len(genome) - (len(genome) % sequence_length)):
        if nucleotide == 0 or nucleotide % sequence_length == 0:
            sequences.append(genome[nucleotide : nucleotide + sequence_length])
    return sequences


def splitted_genome_coordinates(splitted_genome):
    """Coordinates for the genome divided into sequences.
    The function is needed in order to reduce the number of coordinates to construct the z-curve.
    """
    count_x, count_y, count_z = 0, 0, 0
    coordinates = {'x': [], 'y': [], 'z': []}
    for sequence_i, sequence in enumerate(splitted_genome):
        for nucleotide_i in enumerate(sequence):
            if splitted_genome[sequence_i][nucleotide_i[0]] == 'A':
                count_x += 1
                count_y += 1
                count_z += 1
            elif splitted_genome[sequence_i][nucleotide_i[0]] == 'T':
                count_x -= 1
                count_y -= 1
                count_z += 1
            elif splitted_genome[sequence_i][nucleotide_i[0]] == 'G':
                count_x += 1
                count_y -= 1
                count_z -= 1
            elif splitted_genome[sequence_i][nucleotide_i[0]] == 'C':
                count_x -= 1
                count_y += 1
                count_z -= 1
        coordinates['x'].append(count_x * 2)
        coordinates['y'].append(count_y * 2)
        coordinates['z'].append(count_z * 2)
    return coordinates

def main(settings):
    """Analyse input data."""
    input_file = settings["input_files"]
    output_file = settings["output_files"]
    algorithm = settings["algorithm"]

    z_curve_plot = plot_maker(genome_coordinates(delete_description(input_file)))
    # output_z_curve(output_file)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Plots z-curve.")
    parser.add_argument('-i', "--input", help="Genome file in fasta format.", required=True)
    parser.add_argument('-o', "--output", help="File name for the graph")
    parser.add_argument('-f', "--fast", help="Apply a fast algorithm to calculate coordinates", 
                                                                    required=False, default=False)
    args = parser.parse_args()

    input_file = args["input"]
    output_file = args["output"]
    algorithm = args["fast"]

    settings = {
        "input_files": input_file,
        "output_file": output_file,
        "algorithm": algorithm
    }

    main(settings)
