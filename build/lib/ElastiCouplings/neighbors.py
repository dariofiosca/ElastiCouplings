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


def poscar_parser():
    """
    Parses a VASP POSCAR file. This version combines all parsing steps into a single function.
    """
    with open('POSCAR', 'r') as file:
        lines = file.readlines()

    scaling_factor = float(lines[1].strip())
    lattice_vectors = [list(map(float, line.split())) for line in lines[2:5]]
    elements = lines[5].split()
    atom_counts = list(map(int, lines[6].split()))
    coord_type = lines[7].strip()
    is_direct = coord_type.lower() in ['direct', 'd']

    # Parse atomic coordinates
    coordinates_start = 8
    coord = [list(map(float, line.split())) for line in lines[coordinates_start:]]

    return {
        'scaling_factor': scaling_factor,
        'lattice_vectors': lattice_vectors,
        'elements': elements,
        'atom_counts': atom_counts,
        'coordinates': coord,
        'coordinate_type': 'Direct' if is_direct else 'Cartesian'
    }

# def find_neighbors():
#     print('DAJE2')
