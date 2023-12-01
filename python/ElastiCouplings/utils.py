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


def write_q_mesh():
    # Define high symmetry points in reciprocal space
    Gamma = np.array([0, 0, 0])
    K = np.array([3 / 8, 3 / 8, 3 / 4])
    L = np.array([1 / 2, 1 / 2, 1 / 2])
    W = np.array([1 / 2, 1 / 4, 3 / 4])
    X = np.array([1 / 2, 0, 1 / 2])

    # Function to interpolate between two points
    def interpolate(start, end, num_points):
        return np.linspace(start, end, num_points, endpoint=False)

    # Define the path and number of points between each
    num_points = 40  # for example
    path = np.concatenate([
        interpolate(K, Gamma, num_points),
        interpolate(Gamma, L, num_points),
        interpolate(L, W, nuxm_points),
        interpolate(W, X, num_points),
        interpolate(X, K, num_points)
    ])


def load_phonons(filename):
    bands = []
    current_band = []
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            if line.strip() == '':
                if current_band:
                    bands.append(np.array(current_band))
                    current_band = []
            else:
                x, y = map(float, line.split())
                current_band.append([x, y])
    return bands


def print_mat(mat):
    """
    Prints a matrix on terminal

    Inputs
    -------
    Matrix : np.array

    """
    ndim = np.shape(mat)[0]
    for ind in range(ndim):
        str = ''
        for ind1 in range(ndim):
            str += " %6.5f   " % (mat[ind, ind1])
        print(str)
