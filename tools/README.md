This directory contains small python scripts and C++ programs which may be used as subsidiary tools.

## Python scripts 

- displace.py : Generate input structure files of displaced configurations for VASP, Quantum-ESPRESSO, OpenMX, LAMMPS, and xTAPP.
- extract.py : Extract atomic displacements, forces, total energies, and effective charges from output files.
- plotband.py : Visualize phonon bands
- plotdos.py : Visualize phonon DOS and atom-projected phonon DOS
- scph_to_qefc.py : Create a new Quantum-ESPRESSO force constant file (*.fc) with anharmonic correction.

To use the scripts, Python environment (+ Numpy) is necessary. 
Matplotlib is also required for plotband.py and plotdos.py. 

To see available options of each script, please run the script with ``--help`` option.

## C++ scripts (binaries)

- analyze\_phonons.\* : Computes (and prints out) phonon lifetimes, mean-free-path, and (cumulative) thermal-conductivity using the .result file as an input
- qe2alm.\* : Converts a Quantum-ESPRESSO force constant file to the ALAMODE XML format.
- dfc2.\* : Create effective harmonic force constant (ALAMODE XML) from input harmonic force constant (XML format) and PREFIX.scph_dfc2.
- virtual.\* : Performs linear interpolation of force constants

To use these code, please edit the Makefile and do make (or use Cmake).

