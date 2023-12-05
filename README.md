# __ElastiCouplings__

---


This program implements the irreducbile representation projective approach for the calculation of onsite and intersite elastic couplings. The formalism is given in . It is a post-processing program that works on top of a fully converged Density Functional Theory + [phonopy](https://phonopy.github.io/phonopy/index.html) calculation with [VASP](https://www.vasp.at).   


### Installation

You can install this release by using pip

```
pip install .
```
or simply by typing

```
python setup.py install
```


### Usage

Examples are provided in the Example folder

The general workflow is the following:

1. DFPT calculation with VASP with cubic symmetry and global reference frame coinciding with the local octahedral one.  

2. Run phonopy command (for more details have a look at the [phonopy](https://phonopy.github.io/phonopy/index.html) documentation)

```
phonopy --fc vasprun.xml
```

to generate the FORCE\_CONSTANTS matrix, then copy in your folder the POSCAR and FORCE\_CONSTANTS matrix.

3. Now two files are needed for the calculation of the elastic couplings. The first is "nearest\_neighbors.dat" and is a list of atomic sites in the POSCAR file that indicate the position of the corresponding ligand ions separated by a space. Per each line 12 numbers are expected, and are the 6 ligand indices of the centre ion and the 6 of the interacting octahedral centre.  It shall look something like this <br> 
<br> 222 286 158 210 277 150   272 144 304 260 136 294
<br> 222 286 158 210 277 150   269 141 301 257 133 295 <br> 
<br> The second file is "bonds.dat" and contains the connecting lattice vectors that bring the centre octahedra to the interacting one. In our case will look something like <br>
<br>0.5  0.5  0
<br>-0.5 -0.5  0

These files can be both generated in python by using:
```
from  ElastiCouplings import neighbors

neighbors.VASPstruct('O', 'Re', 104)
```

where 'O' is the ligand species, 'Re' is the transition metal one and 104 is the number corresponsing to the centre of the origin-octahedra. 

4. Finally the elastic couplings can be calculated. To do so one has to call in the same folder 

```
from  ElastiCouplings import couplings

couplings.calc_couplings(2)
```

where 62 is the maximum number of nearest-neighbors, next-nearest-neighbors and so on to consider. Bear in mind that the maximum number cannot exceed the number of lines/octahedra considered in you "nearest\_neighbors.dat" and "bonds.dat" files.  

---

author: D. Fiore Mosca

