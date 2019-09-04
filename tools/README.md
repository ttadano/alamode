This directory contains small python scripts and C++ programs which may be used as subsidiary tools.

Each Python scripts does the followings:

* displace.py : script to generate input files of displaced configurations for VASP, Quantum-ESPRESSO, and xTAPP.
* extract.py : script to extract atomic displacements, forces, and total energies from output files.
* plotband.py : script for visualizing phonon bands
* plotdos.py : script for visualizing phonon DOS
* scph_to_qefc.py : script to create a new Quantum-ESPRESSO force constant file (*.fc) with anharmonic correction.

To use the scripts, Python environment (+ Numpy) is necessary.
Matplotlib is also required for plotband.py and plotdos.py. 

Usage of each script may be found in the header part of the source.

Each c++ code does the followings:

* analyze\_phonons.\* : programs for analyzing phonon lifetimes, mean-free-path, 
and (cumulative) thermal-conductivity using the .result file as an input
* qe2alm.\* : program to convert Quantum-ESPRESSO force constant files to 
ALAMODE xml files.

To use these code, please edit the Makefile and do make.

