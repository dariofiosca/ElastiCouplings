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
import itertools


class VASPstruct:
    """
    The VASPstruct class initializes the structure file as in ...

    """

    def __init__(self, ligand, ion, oct):
        self.tol = 2
        self.atom_positions = {}
        self.ligand = ligand
        self.q_types = []
        self.el_names = []

        # New attributes
        self.ntypes = 0
        self.nq = 0
        self.nions = []
        self.type_of_ion = []

        # Attributes for lattice vectors and angles
        self.lattice_vecs = None
        # Bravais matrix
        self.bravais = None
        self.a_brav = None
        self.orth = True
        self.lattype = None
        self.verbosity = 1

        self.poscar_parser()
        self.extend_ligands()
        self.find_nearest_neighbors(oct, ion)

    def poscar_parser(self):
        """
        Parses the input file to extract lattice information, atom positions, and other
        relevant structural data. It updates the class attributes accordingly.

        Attributes Modified
        -------------------
        atom_positions : dict
            Dictionary mapping atom names to their positions.
        el_names : list
            List of element names found in the file.
        q_types : list
            List of atomic positions for each element type.
        nions : list
            List of the number of ions for each element type.
        lattice_vecs : ndarray
            Lattice vectors.
        bravais : ndarray
            The Bravais matrix for the structure.
        ntypes : int
            The number of unique atom types.
        nq : int
            Total number of atoms.
        type_of_ion : list
            Type of ion corresponding to each atom in the system.
        """

        verbosity = 1

        with open('POSCAR', 'r') as file:
            lines = file.readlines()

        scaling_factor = float(lines[1].strip())

        self.bravais = np.array([list(map(float, line.split())) for line in lines[2:5]])
        self.a_brav = self.bravais * scaling_factor
        self.el_names = lines[5].split()
        self.ntypes = list(map(int, lines[6].split()))
        self.nq = sum(self.ntypes)
        self.nions = [i for i in range(1, self.nq + 1)]

        coord_type = lines[7].strip()
        is_direct = coord_type.lower() in ['direct', 'd']

        # Parse atomic coordinates
        coordinates_start = 8
        self.q_types = [np.array(list(map(float, line.split()))) for line in lines[coordinates_start:]]

        if is_direct:
            for i, vec in enumerate(self.q_types):
                self.q_types[i] = (self.bravais.T @ vec.T).T

        self.type_of_ion = [name for count, name in zip(self.ntypes, self.el_names) for _ in range(count)]

        if self.verbosity > 0:
            print("Q-types:\n", self.q_types)
            print("Element names:\n", self.el_names)
            print("Number of types:\n", self.ntypes)
            print("Total number of atoms:\n", self.nq)
            print("Number of each atom:\n", self.nions)
            print("Type of ion:\n", self.type_of_ion)
            print("Bravais Matrix:\n", self.a_brav)

    def extend_ligands(self):

        shift = list(itertools.product(range(-2, 3), repeat=3))
        for i in range(len(shift)):
            for ind, lig in enumerate(self.type_of_ion):
                if lig == self.ligand:
                    T_vec_ini = np.dot(self.a_brav.T, np.array(shift[i]).T).T
                    self.q_types.append(self.q_types[ind] + T_vec_ini)
                    self.nions.append(ind + 1)
                    for k, name in enumerate(self.el_names):
                        if name == self.ligand:
                            self.ntypes[k] += 1

        self.nq = sum(self.ntypes)
        self.type_of_ion = [name for count, name in zip(self.ntypes, self.el_names) for _ in range(count)]

        if self.verbosity > 1:
            print("Q-types:\n", self.q_types)
            print("Element names:\n", self.el_names)
            print("Number of types:\n", self.ntypes)
            print("Total number of atoms:\n", self.nq)
            print("Number of each atom:\n", self.nions)
            print("Type of ion:\n", self.type_of_ion)
            print("Bravais Matrix:\n", self.a_brav)

    def find_nearest_neighbors(self, origin, name):

        r_centre = self.q_types[origin - 1]
        print("ORIGIN ATOM:\n", r_centre)

        distance = []
        vectors = []
        r_index = []

        atom_indices = []

        shift = list(itertools.product(range(-1, 2), repeat=3))
        for i in range(len(shift)):
            for ind, atom in enumerate(self.type_of_ion):
                if atom == name:
                    T_vec_ini = np.dot(self.a_brav.T, np.array(shift[i]).T).T
                    atom_indices.append((self.q_types[ind] + T_vec_ini, ind + 1))

        for ion in atom_indices:
            if ion[1] != (origin):
                R_vec = ion[0] - r_centre  # + ion[0]
                vectors.append(R_vec)
                distance.append(np.linalg.norm(R_vec))
                r_index.append(ion[1])

        R_bas = sorted(zip(distance, vectors, r_index), key=lambda x: x[0])

        origin_octahedron, vector_octahedron = self.find_octahedra(origin, find_ref=True)

        with open('nearest_neighbors.dat', 'w') as file:
            with open('bonds.dat', 'w') as file2:
                for r in R_bas:
                    if np.linalg.norm(r[1]) < 20:
                        strb = ' '.join(map(str, np.round(np.linalg.inv(self.a_brav.T) @ r[1], 4))) + '\n'
                        new_octahedra, new_vec_octahedra = self.find_octahedra(r[2])

                        str1 = ' '.join(map(str, origin_octahedron))
                        str2 = ' '.join(map(str, new_octahedra))

                        # Concatenate the string representations
                        out_str = str1 + '    ' + str2 + '\n'
                        file2.write(strb)
                        file.write(out_str)

    def find_octahedra(self, oct, find_ref=True):

        r_centre = self.q_types[oct - 1]
        octahedra = []
        distance = []
        vectors = []

        for ind, lig in enumerate(self.type_of_ion):
            if lig == self.ligand:
                octahedra.append(self.nions[ind])
                vectors.append(self.q_types[ind] - r_centre)
                distance.append(np.linalg.norm(self.q_types[ind] - r_centre))

        oct_bas = sorted(zip(distance, octahedra, vectors), key=lambda x: x[0])

        seen = set()
        oct_bas_no_dupl = []
        for item in oct_bas:
            # Convert numpy.ndarray to a tuple for hashability
            hash_it = (item[0], item[1], tuple(item[2]))

            if hash_it not in seen:
                seen.add(hash_it)
                oct_bas_no_dupl.append(item)

        octa = []
        for ind in oct_bas_no_dupl:
            if np.round(ind[0] - oct_bas_no_dupl[0][0], self.tol) <= 0.1:
                octa.append(ind)

        if find_ref == True:
            new_octa, vec_octa = self.find_reference_frame(octa)
        else:
            new_octa = []
            vec_octa = []
            for ind in octa:
                new_octa.append(ind[1])
                vec_octa.append(ind[2])

        return new_octa, vec_octa

    def parallel(self, vec1, vec2):

        vec1 = np.around(vec1, decimals=6)
        vec2 = np.around(vec2, decimals=6)

        # Check for zero vectors
        if np.all(vec1 == 0) or np.all(vec2 == 0):
            return False

        # Handle division by zero and check for parallelism
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = vec2 / vec1
            if np.isinf(ratio).any() or np.isnan(ratio).any():
                ratio = np.where(np.isinf(ratio) | np.isnan(ratio), 0, ratio)

        # Check if all non-zero elements are the same and positive
        non_zero_ratios = ratio[ratio != 0]
        return len(np.unique(non_zero_ratios)) == 1 and non_zero_ratios[0] > 0

    def find_reference_frame(self, octa):

        # DEFINITION !!!
        # z -> x -> y -> -z -> -x -> -y
        basis = [np.array([0, 0, 1]), np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, -1]),
                 np.array([-1, 0, 0]), np.array([0, -1, 0])]

        new_octa = []
        vec_octa = []
        logics = [False, False, False, False, False, False]

        for i in range(6):
            for ind in octa:
                vec = ind[2] / np.linalg.norm(ind[2])
                if self.parallel(vec, basis[i]) and logics[i] is False:
                    new_octa.append(ind[1])
                    vec_octa.append(ind[2])
                    logics[i] = True

        return new_octa, vec_octa
