#!/usr/bin/env python
#
# extract.py
#
# Simple script to extract atomic displacements, atomic forces, and
# energies from output files.
# Currently, VASP, Quantum-ESPRESSO, and xTAPP are supported.
#
# Copyright (c) 2014, 2015, 2016 Terumasa Tadano
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#

"""
This python script extracts atomic displacements, atomics forces,
and energies.
"""

from __future__ import print_function
import numpy as np
import optparse
import interface.VASP as vasp
import interface.QE as qe
import interface.xTAPP as xtapp
import interface.OpenMX as openmx
import interface.LAMMPS as lammps

usage = "usage: %prog [options] vasprun*.xml (or *.pw.out or *.str or *.md or *.out)"
parser = optparse.OptionParser(usage=usage)

parser.add_option('--VASP',
                  metavar='orig.POSCAR',
                  help="VASP POSCAR file with equilibrium atomic \
                        positions (default: None)")

parser.add_option('--QE',
                  metavar='orig.pw.in',
                  help="Quantum-ESPRESSO input file with equilibrium\
                  atomic positions (default: None)")

parser.add_option('--xTAPP',
                  metavar='orig.cg',
                  help="xTAPP CG file with equilibrium atomic \
                        positions (default: None)")

parser.add_option('--LAMMPS',
                  metavar='orig.lammps',
                  help="LAMMPS structure file with equilibrium atomic positions (default: None)")

parser.add_option('--OpenMX',
                  metavar='orig.dat',
                  help="OpenMX dat file with equilibrium atomic \
                        positions (default: None)")

parser.add_option('--get',
                  default="disp-force",
                  help="specify which quantity to extract. \
                        Available options are 'disp-force', 'disp', 'force' and 'energy'.")

parser.add_option('--unit',
                  action="store",
                  type="string",
                  dest="unitname",
                  default="Rydberg",
                  help="print atomic displacements and forces in units of UNIT. \
                        Available options are 'eV', 'Rydberg' (default), and 'Hartree'.")

parser.add_option('--offset',
                  help="Specify an output file (either *.xml, *.pw.out, or *.str) of an\
 equilibrium structure to subtract residual forces, displacements, or energies.")

parser.add_option('--emin',
                  default=None,
                  type="float",
                  help="Lower bound of the energy filter (eV) used for selecting output structures.\
                        Available only in the VASP parser.")

parser.add_option('--emax',
                  default=None,
                  type="float",
                  help="Upper bound of the energy filter (eV) used for selecting output structures.\
                        Available only in the VASP parser.")

# Main

if __name__ == "__main__":

    options, args = parser.parse_args()
    file_results = args[0:]

    if len(file_results) == 0:
        print("Usage: extract.py [options] vasprun*.xml \
(or *.pw.out, or *.str)")
        print()
        print("For details of available options, please type\n\
$ python displace.py -h")
        exit(1)

    # Check the calculator option

    conditions = [options.VASP is None,
                  options.QE is None,
                  options.xTAPP is None,
                  options.LAMMPS is None,
                  options.OpenMX is None]

    if conditions.count(True) == len(conditions):
        print(
            "Error : Either --VASP, --QE, --xTAPP, --LAMMPS, --OpenMX option must be given.")
        exit(1)

    elif len(conditions) - conditions.count(True) > 1:
        print("Error : --VASP, --QE, --xTAPP, --LAMMPS, and --OpenMX cannot be given simultaneously.")
        exit(1)

    elif options.VASP:
        code = "VASP"
        file_original = options.VASP

    elif options.QE:
        code = "QE"
        file_original = options.QE

    elif options.xTAPP:
        code = "xTAPP"
        file_original = options.xTAPP

    elif options.LAMMPS:
        code = "LAMMPS"
        file_original = options.LAMMPS

    elif options.OpenMX:
        code = "OpenMX"
        file_original = options.OpenMX

    # Check output option
    str_get = options.get.lower()
    if str_get not in ["disp-force", "disp", "force", "energy"]:
        print("Error: Please specify which quantity to extract by the --get option.")
        exit(1)

    print_disp = False
    print_force = False
    print_energy = False

    if str_get == "disp-force":
        print_disp = True
        print_force = True
    elif str_get == "disp":
        print_disp = True
    elif str_get == "force":
        print_force = True
    elif str_get == "energy":
        print_energy = True

    # Check unit option
    str_unit = options.unitname.lower()
    if str_unit in ["ev", "electron_volt"]:
        str_unit = "ev"
    elif str_unit in ["ry", "ryd", "rydberg"]:
        str_unit = "rydberg"
    elif str_unit in ["ha", "hartree"]:
        str_unit = "hartree"
    else:
        print("Error: Invalid unit name : %s" % options.unitname)

    # Print data
    if code == "VASP":
        vasp.parse(file_original, file_results,
                   options.offset, str_unit,
                   print_disp, print_force, print_energy,
                   options.emin, options.emax)

    elif code == "QE":
        qe.parse(file_original, file_results,
                 options.offset, str_unit,
                 print_disp, print_force, print_energy)

    elif code == "xTAPP":
        xtapp.parse(file_original, file_results,
                    options.offset, str_unit,
                    print_disp, print_force, print_energy)

    elif code == "OpenMX":
        openmx.parse(file_original, file_results,
                     options.offset, str_unit,
                     print_disp, print_force, print_energy)

    elif code == "LAMMPS":
        lammps.parse(file_original, file_results,
                     options.offset, str_unit,
                     print_disp, print_force, print_energy)
