This directory contains small python scripts and C++ programs which may be used as subsidiary tools.

Each Python scripts does the followings:

* displace.py : script to generate input files of displaced configurations for VASP, Quantum-ESPRESSO, and xTAPP.
* extract.py : script to extract atomic displacements, forces, and total energies from output files.
* plotband.py : script for visualizing phonon bands
* plotdos.py : script for visualizing phonon DOS
* analyze\_phonons.\* : programs for analyzing phonon lifetimes, mean-free-path, and (cumulative) thermal-conductivity 
using the .result file as an input

To use the scripts, Python environment (+ Numpy) is necessary.
Matplotlib is also required for plotband.py and plotdos.py. 

Usage of each script may be found in the header part of the source.


The following scripts are deprecated, and will be removed from the package in a future release.
Please use displace.py or extract.py.

* genPOSCAR.py : script for VASP users which makes POSCAR files of displaced configurations
* genCG.py : script for xTAPP users which makes \*.cg files of displaced configurations
* extVASP.py : script for VASP users which extvacts atomic displacements/forces from vasprun.xml files
* extTAPP.py : script for xTAPP users which extvacts atomic displacements/forces from \*.str files


