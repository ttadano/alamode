#!/usr/bin/env python
#
# displace.py
#
# Simple script to generate input files of given displacement patterns.
# Currently, VASP, Quantum-ESPRESSO, LAMMPS, OpenMX, and xTAPP are supported.
#
# Copyright (c) 2014-2020 Terumasa Tadano
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#

"""
Input structure file generator for displaced configurations.
"""

from __future__ import print_function
import argparse
from interface.VASP import VaspParser
from interface.QE import QEParser
from interface.xTAPP import XtappParser
from interface.OpenMX import OpenmxParser
from interface.LAMMPS import LammpsParser
from GenDisplacement import AlamodeDisplace

parser = argparse.ArgumentParser()

parser.add_argument('--mag',
                    type=float, default=0.02,
                    help="Magnitude of displacement in units of \
                        Angstrom (default: 0.02)")

parser.add_argument('--prefix',
                    type=str, default="disp",
                    help="Prefix of the files to be created. (default: disp)")

parser.add_argument('--QE',
                    metavar='supercell.pw.in',
                    help="Quantum-ESPRESSO input file with equilibrium atomic positions (default: None)")

parser.add_argument('--VASP',
                    metavar='SPOSCAR',
                    help="VASP POSCAR file with equilibrium atomic \
                        positions (default: None)")

parser.add_argument('--xTAPP',
                    metavar='supercell.cg',
                    help="xTAPP CG file with equilibrium atomic \
                        positions (default: None)")

parser.add_argument('--LAMMPS',
                    metavar='supercell.lammps',
                    help="LAMMPS structure file with equilibrium atomic positions (default: None)")

parser.add_argument('--OpenMX',
                    metavar='supercell.dat',
                    help="dat file with equilibrium atomic \
                        positions (default: None)")

parser.add_argument('-pf', '--pattern_file', metavar='prefix.pattern_*', type=str, nargs='+',
                    help="ALM pattern file(s) generated with MODE = suggest")

parser.add_argument('--random', action="store_true", dest="random", default=False,
                    help="Generate randomly-displaced structures.")

parser.add_argument('--temp', type=float, default=None,
                    help="Target temperature of the random distribution of \
                        Q (default: None). Used if --MD is not given.")

parser.add_argument('-e', '--every', type=str, default="50", metavar='start:end:interval',
                    help="Specify the range and interval of data sampling. "
                         "--every=1:1000:10 means sampling one structure for every 10 snapshots"
                         "from the 1st step to the 1000th step.Used if --MD is given. "
                         "(default: 50)")

parser.add_argument('-md', '--load_mddata', type=str, nargs='+',
                    help="Specify the file(s) containing displacements of MD trajectories.")

parser.add_argument('--prim', type=str, default=None,
                    help="Specify the file containing structure data of the primitive lattice.")

parser.add_argument('--evec', type=str, default=None,
                    help="Specify the file containing harmonic eigenvalues and eigenvectors.")

parser.add_argument('-nd', '--num_disp', type=int, default=1,
                    help="Specify the number of displacement patterns.")

parser.add_argument('-cl', '--classical', action="store_true", dest="classical", default=False,
                    help="Use classical expectation value for <Q^2>.")

parser.add_argument('-p', '--print', action="store_true", dest="print_disp", default=False,
                    help="Print displacements to stdout")

parser.add_argument('--pes', type=str, default=None,
                    help="Specify the target mode to compute PES.")

parser.add_argument('--Qrange', type=str, default=None,
                    help='Range of normal coordinate Q in units of eV/Ang.')


def check_options(args):
    conditions = [args.VASP is None,
                  args.QE is None,
                  args.xTAPP is None,
                  args.LAMMPS is None,
                  args.OpenMX is None]

    if conditions.count(True) == len(conditions):
        raise RuntimeError(
            "Error : Either --VASP, --QE, --xTAPP, --LAMMPS, "
            "--OpenMX option must be given.")

    elif len(conditions) - conditions.count(True) > 1:
        raise RuntimeError("Error : --VASP, --QE, --xTAPP, --LAMMPS, and "
                           "--OpenMX cannot be given simultaneously.")

    elif args.VASP:
        code = "VASP"
        struct_format = "VASP POSCAR"
        str_outfiles = "%s{counter}.POSCAR" % args.prefix
        file_original = args.VASP

    elif args.QE:
        code = "QE"
        struct_format = "Quantum-ESPRESSO pw.in format"
        str_outfiles = "%s{counter}.pw.in" % args.prefix
        file_original = args.QE

    elif args.xTAPP:
        code = "xTAPP"
        struct_format = "xTAPP cg format"
        str_outfiles = "%s{counter}.cg" % args.prefix
        file_original = args.xTAPP

    elif args.LAMMPS:
        code = "LAMMPS"
        struct_format = "LAMMPS structure format"
        str_outfiles = "%s{counter}.lammps" % args.prefix
        file_original = args.LAMMPS

    elif args.OpenMX:
        code = "OpenMX"
        struct_format = "OpenMX dat format"
        str_outfiles = "%s{counter}.dat" % args.prefix
        file_original = args.OpenMX

    return code, file_original, struct_format, str_outfiles


def get_code_object(code):
    if code == "VASP":
        return VaspParser()

    if code == "QE":
        return QEParser()

    if code == "OpenMX":
        return OpenmxParser()

    if code == "xTAPP":
        return XtappParser()

    if code == "LAMMPS":
        return LammpsParser()


if __name__ == '__main__':
    print("*****************************************************************")
    print("    displace.py --  Generator of displaced configurations        ")
    print("                      Version. 1.2.0                             ")
    print("*****************************************************************")
    print("")

    args = parser.parse_args()
    file_pattern = args.pattern_file

    code, file_original, struct_format, str_outfiles = check_options(args)

    codeobj = get_code_object(code)
    codeobj.load_initial_structure(file_original)

    print(" Output format                  : %s" % struct_format)
    print(" Structure before displacements : %s" % file_original)
    print(" Output file names              : %s" % str_outfiles)
    print(" Magnitude of displacements     : %s Angstrom" % args.mag)
    print(" Number of atoms                : %i" % codeobj.nat)
    print("")

    dispobj = AlamodeDisplace("fd", codeobj,
                              file_primitive=args.prim,
                              file_evec=args.evec)
    header_list, disp_list = dispobj.generate(file_pattern=file_pattern,
                                              file_mddata=args.load_mddata,
                                              option_every=args.every,
                                              magnitude=args.mag,
                                              number_of_displacements=args.num_disp,
                                              temperature=args.temp)

    codeobj.generate_structures(args.prefix, header_list, disp_list)
    print(" Number of displacements        : %i" % len(disp_list))
    print("-----------------------------------------------------------------")
    print("")
    print("All input files are created.")
