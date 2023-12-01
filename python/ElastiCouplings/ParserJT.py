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


def jahn_teller_ops():
    # Initialize Jahn-Teller

    Q1 = np.zeros((6, 1, 3, 3))

    Q2 = np.zeros((6, 1, 3, 3))
    Q3 = np.zeros((6, 1, 3, 3))

    Qyz = np.zeros((6, 1, 3, 3))
    Qxz = np.zeros((6, 1, 3, 3))
    Qxy = np.zeros((6, 1, 3, 3))

    # ---------> Q1

    sq16 = 1 / np.sqrt(6)
    q1_ind = [sq16, -sq16, sq16, -sq16, sq16, -sq16]
    q1_val = [(1, 0, 0, 0), (4, 0, 0, 0), (2, 0, 1, 1), (5, 0, 1, 1), (0, 0, 2, 2), (3, 0, 2, 2)]
    for i, ind in enumerate(q1_val):
        Q1[ind] = q1_ind[i]

    # ---------> Eg
    # Q2
    q2_ind = [0.5, -0.5, -0.5, 0.5]
    q2_val = [(1, 0, 0, 0), (4, 0, 0, 0), (2, 0, 1, 1), (5, 0, 1, 1)]
    for i, ind in enumerate(q2_val):
        Q2[ind] = q2_ind[i]
    # Q3
    sq3 = 1 / np.sqrt(3)
    q3_ind = [sq3, -sq3, -0.5 * sq3, 0.5 * sq3, -0.5 * sq3, 0.5 * sq3]
    q3_val = [(0, 0, 2, 2), (3, 0, 2, 2), (1, 0, 0, 0), (4, 0, 0, 0), (2, 0, 1, 1), (5, 0, 1, 1)]
    for i, ind in enumerate(q3_val):
        Q3[ind] = q3_ind[i]

    #  ---------> T2g
    # Qyz
    qyz_ind = [0.5, -0.5, 0.5, -0.5]
    qyz_val = [(2, 0, 2, 2), (5, 0, 2, 2), (0, 0, 1, 1), (3, 0, 1, 1)]
    for i, ind in enumerate(qyz_val):
        Qyz[ind] = qyz_ind[i]

    # Qxz
    qxz_val = [(0, 0, 0, 0), (3, 0, 0, 0), (1, 0, 2, 2), (4, 0, 2, 2)]
    for i, ind in enumerate(qxz_val):
        Qxz[ind] = qyz_ind[i]

    # Qxy
    qxy_val = [(1, 0, 1, 1), (4, 0, 1, 1), (2, 0, 0, 0), (5, 0, 0, 0)]
    for i, ind in enumerate(qxy_val):
        Qxy[ind] = qyz_ind[i]

    Q1 = Q1.reshape(18, 3)
    Q2 = Q2.reshape(18, 3)
    Q3 = Q3.reshape(18, 3)
    Qxz = Qxz.reshape(18, 3)
    Qyz = Qyz.reshape(18, 3)
    Qxy = Qxy.reshape(18, 3)

    Q = np.concatenate((Qxy, Qyz, Q3, Qxz, Q2), axis=1)

    return Q
