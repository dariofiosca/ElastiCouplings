# ------------------------------------------------------------------------------------#
# ElastiCouplings library
# Written by  : Dario Fiore Mosca (Ecole Polytechnique) 2023
# Email: dario.fiore.mosca@univie.ac.at
# ------------------------------------------------------------------------------------#
#
#    ElastiCouplings implements the Elastic Coupling calculation of Ref.
#
# ------------------------------------------------------------------------------------#
from ElastiCouplings.utils import *
import matplotlib.pyplot as plt
import numpy as np

# Step 1: Parse the bands_new.dat file
x_axis_labels = []
irreps_labels = []

with open("irreps_labels", "r") as file:
    for line in file:
        if line.strip():
            x, y, label = line.split()
            x = float(x)
            y = float(y)
            irreps_labels.append((x, y, label))

plt.figure(figsize=(6, 6))

filename_phonopy = 'bands_new_2.dat'
bands_phonopy = load_phonons(filename_phonopy)
# Step 2: Plot the data
for index, band in enumerate(bands_phonopy):
    plt.plot(band[:, 0], band[:, 1], linewidth=3.1, color='black')
    plt.plot(band[:, 0], band[:, 1], linewidth=3, color='grey')
    if index == 0:
        plt.plot(band[:, 0], band[:, 1], linewidth=3, color='grey', label='Phonopy')

filename_projected = 'phonons_nn.dat'
bands_projected = load_phonons(filename_projected)
for index, band in enumerate(bands_projected):
    plt.plot(bands[0][:, 0], band[:, 1], linewidth=3.1, color='black')
    plt.plot(bands[0][:, 0], band[:, 1], linewidth=3, color='darkblue')
    if index == 0:
        plt.plot(bands[0][:, 0], band[:, 1], linewidth=3, color='darkblue', label='Onsite + NN')

# Plotting the labels from irreps_labels
for x, y, label in irreps_labels:
    plt.text(x, y, label, fontsize=14, ha='right', va='bottom')

# Set custom ticks and labels on the x-axis
tick_positions = [0.00000000, 0.07710320, 0.17990740, 0.26384670, 0.32320080]
tick_labels = ['K', 'G', 'L', 'W', 'X']
plt.xticks(tick_positions, tick_labels)

plt.tick_params(axis='both', which='major', labelsize=14)
plt.xlim(0.0, 0.32320080)
plt.ylim(8, 22)
plt.xlabel('Wave Vector', fontsize=14)
plt.ylabel('Frequency (THz)', fontsize=14)
plt.title('Phonon Band Structure')
plt.grid(True, axis='x')
plt.legend()

plt.savefig('phonons_bcro_test.pdf', dpi=1000)

plt.show()
