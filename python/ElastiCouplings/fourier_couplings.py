# ------------------------------------------------------------------------------------#
# ElastiCouplings library
# Written by  : Dario Fiore Mosca (Ecole Polytechnique) 2023
# Email: dario.fiore.mosca@univie.ac.at
# ------------------------------------------------------------------------------------#
#
#    ElastiCouplings implements the Elastic Coupling calculation of Ref.
#
# ------------------------------------------------------------------------------------#
import numpy as np
import sys

# Initialization
f_meV = 16.159  # convertion eigenvalues to meV for oxygen mass
f_cm1 = 130.34  # convertion eigenvalues to cm-1 for oxygen mass
filename = sys.argv[1]

q_vectors = []
M_prim = np.array([[0., 0.5, 0.5], [0.5, 0., 0.5], [0.5, 0.5, 0.]])

# Open the file and read each line
with open("q_mesh", "r") as file:
    for line in file:
        vector_components = [float(x) for x in line.split()]
        q_vectors.append(np.array(vector_components))

R_values = []  # List to store R values
matrices = []  # List to store 5x5 matrices

# Read File
with open(filename, 'r') as file:
    lines = file.readlines()

i = 0
while i < len(lines):
    line = lines[i].strip()
    if line.startswith('R ='):
        # Extract and clean R values
        R = line[4:].strip('[] ').split()
        R = [float(val) for val in R]
        # Attention here the matrix is projected in the primitive basis
        R_values.append(M_prim @ R)

        i += 1
        matrix = []

        # Read the next 5 lines to fill the matrix
        for _ in range(5):
            matrix_line = lines[i].strip().split()
            matrix.append([complex(x.replace('i', 'j')) for x in matrix_line])  # Convert to complex
            i += 1

        matrices.append(np.array(matrix, dtype=complex))

    else:
        i += 1

# TODO Read the t2g_onsite and eg_onsite from Elastic_constants.dat
t2g_on_site = 8.23989484147813  # on-site in eV/Ang^2
eg_on_site = 20.566994115031463
on_site = np.zeros((5, 5), dtype=complex)

# on_site[0,0]=q1_on_site
for i in [0, 1, 3]:
    on_site[i, i] = t2g_on_site
for i in [2, 4]:
    on_site[i, i] = eg_on_site

for daje in range(5):
    print(' ')
    for q_ind, q in enumerate(q_vectors):
        H = np.zeros((5, 5), dtype=complex)
        H += on_site
        for r_ind, r in enumerate(R_values):
            # print(np.linalg.inv(M_prim) @ r)
            phase = 2 * np.pi * np.dot(q, r)
            exp_fac = np.exp(complex(0.0, phase))
            H += exp_fac * matrices[r_ind]
        e, v = np.linalg.eigh(H)
        print("%s %10.3f" % (q_ind, np.sqrt(e[daje]) * f_cm1 * 0.030))
