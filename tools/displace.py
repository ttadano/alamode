#!/usr/bin/env python
#
# displace.py
#
# Simple script to generate input files of given displacement patterns.
# Currently, VASP, Quantum-ESPRESSO, and xTAPP are supported.
#
# Copyright (c) 2014 Terumasa Tadano
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#

"""
Input file generator for displaced configurations.
"""

from __future__ import print_function
import argparse
import numpy as np
import interface.VASP as vasp
import interface.QE as qe
import interface.xTAPP as xtapp
import interface.OpenMX as openmx
import interface.LAMMPS as lammps

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

parser.add_argument('pattern_file', metavar='prefix.pattern_*', type=str, nargs='+',
                    help="ALM pattern file(s) generated with MODE = suggest")



def parse_displacement_patterns(files_in):

    pattern = []

    for file in files_in:
        pattern_tmp = []

        f = open(file, 'r')
        tmp, basis = f.readline().rstrip().split(':')
        if basis == 'F':
            print("Warning: DBASIS must be 'C'")
            exit(1)

        while True:
            line = f.readline()

            if not line:
                break

            line_split_by_colon = line.rstrip().split(':')
            is_entry = len(line_split_by_colon) == 2

            if is_entry:
                pattern_set = []
                natom_move = int(line_split_by_colon[1])
                for i in range(natom_move):
                    disp = []
                    line = f.readline()
                    line_split = line.rstrip().split()
                    disp.append(int(line_split[0]))
                    for j in range(3):
                        disp.append(float(line_split[j + 1]))

                    pattern_set.append(disp)
                pattern_tmp.append(pattern_set)

        print("File %s containts %i displacement patterns"
              % (file, len(pattern_tmp)))

        for entry in pattern_tmp:
            if entry not in pattern:
                pattern.append(entry)

        f.close()

    print("Number of unique displacement patterns = %d" % len(pattern))

    return pattern


def char_xyz(entry):

    if entry % 3 == 0:
        return 'x'
    if entry % 3 == 1:
        return 'y'
    if entry % 3 == 2:
        return 'z'


def gen_displacement(counter_in, pattern, disp_mag, nat, invlavec):

    poscar_header = "Disp. Num. %i" % counter_in
    poscar_header += " ( %f Angstrom" % disp_mag

    disp = np.zeros((nat, 3))

    for displace in pattern:
        atom = displace[0] - 1

        poscar_header += ", %i : " % displace[0]

        str_direction = ""

        for i in range(3):
            if abs(displace[i + 1]) > 1.0e-10:
                if displace[i + 1] > 0.0:
                    str_direction += "+" + char_xyz(i)
                else:
                    str_direction += "-" + char_xyz(i)

            disp[atom][i] += displace[i + 1] * disp_mag

        poscar_header += str_direction

    poscar_header += ")"

    if invlavec is not None:
        for i in range(nat):
            disp[i] = np.dot(disp[i], invlavec.T)

    return poscar_header, disp


def get_number_of_zerofill(npattern):

    nzero = 1

    while True:
        npattern //= 10

        if npattern == 0:
            break

        nzero += 1

    return nzero


def check_options(args):

    conditions = [args.VASP is None,
                  args.QE is None,
                  args.xTAPP is None,
                  args.LAMMPS is None,
                  args.OpenMX is None]

    if conditions.count(True) == len(conditions):
        print(
            "Error : Either --VASP, --QE, --xTAPP, --LAMMPS, --OpenMX option must be given.")
        exit(1)

    elif len(conditions) - conditions.count(True) > 1:
        print("Error : --VASP, --QE, --xTAPP, --LAMMPS, and --OpenMX cannot be given simultaneously.")
        exit(1)

    elif args.VASP:
        code = "VASP"
        struct_format = "VASP POSCAR"
        str_outfiles = "%s{counter}.POSCAR" % args.prefix
        file_original = args.VASP
        suffix = None

    elif args.QE:
        code = "QE"
        struct_format = "Quantum-ESPRESSO pw.in format"
        str_outfiles = "%s{counter}.pw.in" % args.prefix
        file_original = args.QE
        suffix = "pw.in"

    elif args.xTAPP:
        code = "xTAPP"
        struct_format = "xTAPP cg format"
        str_outfiles = "%s{counter}.cg" % args.prefix
        file_original = args.xTAPP
        suffix = None

    elif args.LAMMPS:
        code = "LAMMPS"
        struct_format = "LAMMPS structure format"
        str_outfiles = "%s{counter}.lammps" % args.prefix
        file_original = args.LAMMPS
        suffix = None

    elif args.OpenMX:
        code = "OpenMX"
        struct_format = "OpenMX dat format"
        str_outfiles = "%s{counter}.dat" % args.prefix
        file_original = args.OpenMX
        suffix = None

    return code, file_original, struct_format, str_outfiles, suffix


if __name__ == '__main__':

    print("*****************************************************************")
    print("    displace.py --  Generator of displaced cofigurations         ")
    print("*****************************************************************")
    print("")

    args = parser.parse_args()
    file_pattern = args.pattern_file

    code, file_original, struct_format, str_outfiles, suffix = check_options(args)
    disp_length = args.mag
    prefix = args.prefix

    # Read the original file
    if code == "VASP":
        aa, aa_inv, elems, nats, x_frac = vasp.read_POSCAR(file_original)
        nat = np.sum(nats)

    elif code == "QE":
        list_namelist, list_ATOMIC_SPECIES, \
            list_K_POINTS, list_CELL_PARAMETERS, list_OCCUPATIONS, \
            nat, lavec, kd_symbol, x_frac, aa_inv = qe.read_original_QE(
                file_original)

    elif code == "xTAPP":
        str_header, nat, nkd, aa, aa_inv, x_frac, kd \
            = xtapp.read_CG(file_original)
        suffix = "cg"

    elif code == "LAMMPS":
        common_settings, nat, x_cart, kd, charge \
            = lammps.read_lammps_structure(file_original)
        aa_inv = None

    elif code == "OpenMX":
        aa, aa_inv, nat, x_frac = openmx.read_OpenMX_input(file_original)

    print(" Output format                  : %s" % struct_format)
    print(" Structure before displacements : %s" % file_original)
    print(" Output file names              : %s" % str_outfiles)
    print(" Magnitude of displacements     : %s Angstrom" % disp_length)
    print(" Number of atoms                : %i" % nat)
    print("")
    print("-----------------------------------------------------------------")

    disp_pattern = parse_displacement_patterns(file_pattern)
    nzerofills = get_number_of_zerofill(len(disp_pattern))
    counter = 0

    for pattern in disp_pattern:
        counter += 1
        header, disp = gen_displacement(counter, pattern, disp_length,
                                        nat, aa_inv)

        if code == "VASP":
            vasp.write_POSCAR(prefix, counter, header, nzerofills,
                              aa, elems, nats, disp, x_frac)

        elif code == "QE":
            qe.generate_QE_input(prefix, suffix, counter, nzerofills, list_namelist,
                                 list_ATOMIC_SPECIES, list_K_POINTS,
                                 list_CELL_PARAMETERS, list_OCCUPATIONS,
                                 nat, kd_symbol, x_frac, disp)

        elif code == "xTAPP":
            nsym = 1
            symop = []
            symop.append([1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0])
            denom_tran = 1
            has_inv = 0

            xtapp.gen_CG(prefix, suffix, counter, nzerofills, str_header, nat, kd,
                         x_frac, disp, nsym, symop, denom_tran, has_inv)

        elif code == "LAMMPS":
            lammps.write_lammps_structure(prefix, counter, header, nzerofills,
                                          common_settings, nat, kd, x_cart, disp, charge)

        elif code == "OpenMX":
            openmx.write_OpenMX_input(
                prefix, counter,  nzerofills, disp, aa, file_original)

    print("")
    print("All input files are created.")
