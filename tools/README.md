This directory contains small python scripts and C++ programs which may be used as subsidiary tools.

Each Python scripts does the followings:

* genPOSCAR.py : script for VASP users which makes POSCAR files of displaced configurations
* extVASP.py : script for VASP users which extvacts atomic displacements/forces from vasprun.xml files
* plotband.py : script for visualizing phonon bands
* plotdos.py : script for visualizing phonon DOS
* analyze_phonons.* : programs for analyzing phonon lifetimes, mean-free-path, and thermal-conductivity 
using the .result file as an input

To use the scripts, Python environment is necessarily including Numpy.
To use plotband.py or plotdos.py, Matplotlib is also required. 

Usage of each script may be found in the header part of the source.