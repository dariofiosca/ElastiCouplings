# __ElastiCouplings__

---


This program implements the irreducbile representation projective approach for the calculation of onsite and intersite elastic couplings. The formalism is given in . It is a post-processing program that works on top of a fully converged Density Functional Theory + [phonopy](https://phonopy.github.io/phonopy/formulation.html) calculation with [VASP](https://www.vasp.at).   


### Installation

You can install this release by using pip

```
pip install .
```

### Usage

Examples for how to use it are provided in the Example folder

The general workflow is the following:

1. DFPT calculation with VASP and large unit cell

2. Using phonopy commands (for more details have a look at the phonopy documentation)

```
phonopy --fc vasprun.xml
```
for creating the FORCE\_CONSTANTS matrix.

3. Copy in your folder the POSCAR and FORCE constant matrix. 

...


---

author: D. Fiore Mosca

