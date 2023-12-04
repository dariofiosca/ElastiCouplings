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
from ElastiCouplings.parser import *
from ElastiCouplings.utils import *

# Define the file path
file_path = 'FORCE_CONSTANTS'


# Define the indices list
def calc_couplings(NN):
    # Initialise the nearest neighbor indices between atom 1 (centre) and 2 (NN)
    arr = np.loadtxt('nearest_neighbors.dat')
    at_12 = np.zeros((NN, 12), dtype=float)
    for at in range(NN):
        at_12[at, :] = arr[at, :]

    # Initialise the different bonds
    arr = np.loadtxt('bonds.dat')
    bonds = np.zeros((NN, 3), dtype=float)
    for at in range(NN):
        bonds[at, :] = arr[at, :]

    for at in range(NN):
        # Initialize the smaller matrix
        s_mat = np.zeros((12, 12, 3, 3))

        # Initialize the row and column indices
        row_ind = None
        col_ind = None

        # Read the file line by line
        with open(file_path, 'r') as file:
            # Skip the first line with total rows and columns information
            next(file)

            # Iterate through the remaining lines
            for line in file:
                # Split the line into tokens
                tokens = line.split()

                # Check if the line contains the index of the row or column
                if len(tokens) == 2:
                    row_ind, col_ind = int(tokens[0]), int(tokens[1])

                # Check if the row or column index is in the at_1 list
                if row_ind in at_12[at, :].tolist() or col_ind in at_12[at, :].tolist():
                    # Read the next three lines to populate the submatrix
                    submatrix = [[float(value) for value in next(file).split()] for _ in range(3)]

                    # Find the index of the row and column in the at_1 list
                    row_at_1_index = at_12[at, :].tolist().index(row_ind) if row_ind in at_12[at, :].tolist() else None
                    col_at_1_index = at_12[at, :].tolist().index(col_ind) if col_ind in at_12[at, :].tolist() else None

                    # Assign the submatrix to the corresponding position in the smaller matrix
                    if row_at_1_index is not None and col_at_1_index is not None:
                        s_mat[row_at_1_index, col_at_1_index] = submatrix
                        # print(row_ind, col_ind)
                        # print(row_at_1_index, col_at_1_index)
                        # print_arr(np.array(submatrix))

        # Print the smaller matrix
        original_shape = s_mat.shape
        new_shape = (original_shape[0] * original_shape[2], original_shape[1] * original_shape[3])

        # Reshape and transpose the matrix
        s_mat = s_mat.transpose(0, 2, 1, 3).reshape(new_shape)

        # Define Jahn-Teller operators
        Q = jahn_teller_ops()

        # On-site 1
        zeros_array = np.zeros((18, 3))

        names = ['Qxy', 'Qyz', 'Q3', 'Qxz', 'Q2']

        jt_at_1 = []
        jt_at_2 = []
        jt_at_12 = np.zeros((5, 5))

        for i, val in enumerate(names):
            mat = np.concatenate((Q[:, i * 3:(i * 3 + 3)], zeros_array), axis=0)
            nmat = np.array(mat[:, 0] + mat[:, 1] + mat[:, 2])
            jt_1 = nmat.transpose() @ s_mat @ nmat
            jt_at_1.append(jt_1)

        for i, val in enumerate(names):
            mat = np.concatenate((zeros_array, Q[:, i * 3:(i * 3 + 3)]), axis=0)
            nmat = np.array(mat[:, 0] + mat[:, 1] + mat[:, 2])
            jt_2 = nmat.transpose() @ s_mat @ nmat
            jt_at_2.append(jt_2)

        if at == 0:
            print('\nOn-site 1')
            for i, val in enumerate(names):
                print(val, ' : ', jt_at_1[i])

            print('\nOn-site 2')
            for i, val in enumerate(names):
                print(val, ' : ', jt_at_2[i])
            print('\nTwo Sites ')

        for i, val1 in enumerate(names):
            for j, val2 in enumerate(names):
                mat = np.concatenate((Q[:, i * 3:(i * 3 + 3)], Q[:, j * 3:(j * 3 + 3)]), axis=0)
                nmat = np.array(mat[:, 0] + mat[:, 1] + mat[:, 2])
                jt_12 = nmat.transpose() @ s_mat @ nmat
                jt_at_12[i, j] = jt_12 - jt_at_1[i] - jt_at_2[j]

        print('\nR = ', bonds[at, :])
        print_mat(np.round(jt_at_12, 5))
