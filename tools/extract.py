#!/usr/bin/env python
#
# extract.py
#
# Simple script to extract atomic displacements, atomic forces, and
# energies from output files.
# Currently, VASP, Quantum-ESPRESSO, and xTAPP are supported.
#
# Copyright (c) 2014 Terumasa Tadano
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
"""
This python script extracts atomic displacements, atomics forces,
and energies.
"""

import numpy as np
import optparse
from displace import *

try:
    try:
        # cElementTree on Python 2.5+
        import xml.etree.cElementTree as etree
    except ImportError:
        # ElementTree on Python 2.5+
        import xml.etree.ElementTree as etree
except ImportError:
    try:
        # cElementTree
        import cElementTree as etree
    except ImportError:
        # ElementTree
        import elementtree.ElementTree as etree


usage = "usage: %prog [options] vasprun*.xml (or *.pw.out or *.str)"
parser = optparse.OptionParser(usage=usage)

parser.add_option('--unit',
                  action="store",
                  type="string",
                  dest="unitname",
                  default="Rydberg",
                  help="print atomic displacements and forces in units of UNIT. \
                        Available options are eV, Rydberg (default), and Hartree.")

parser.add_option('--get',
                  help="specify which quantity to extract. \
                        Available options are disp, force and energy.")

parser.add_option('--QE',
                  metavar='orig.pw.in',
                  help="Quantum-ESPRESSO input file with equilibrium\
 atomic positions (default: None)")
parser.add_option('--VASP',
                  metavar='orig.POSCAR',
                  help="VASP POSCAR file with equilibrium atomic \
                        positions (default: None)")
parser.add_option('--xTAPP',
                  metavar='orig.cg',
                  help="xTAPP CG file with equilibrium atomic \
                        positions (default: None)")

# Functions for VASP

def print_displacements_VASP(xml_files, lavec, nat, x0, 
                             require_conversion, conversion_factor):

    x0 = np.round(x0, 8)
    x = np.zeros((nat, 3))

    lavec_transpose = lavec.transpose()

    for search_target in xml_files:

        xml = etree.parse(search_target)
        root = xml.getroot()

        for elems in root.findall('calculation/structure/varray'):
            str_coord = []
            for elems2 in elems.findall('v'):
                str_coord.append(elems2.text)

            n = len(str_coord)

            for i in range(n):
                x[i, :] = [float(t) for t in str_coord[i].split()]

            x = np.matrix(x)
            disp = x - x0

            for i in range(n):
                disp[i, :] = [refold(disp[i, j]) for j in range(3)]

            disp = disp * lavec_transpose

            if require_conversion:
                for i in range(n):
                    disp[i, :] *= conversion_factor

            for i in range(n):
                print "%15.7F %15.7F %15.7F" % (disp[i, 0],
                                                disp[i, 1],
                                                disp[i, 2])


def print_atomicforces_VASP(xml_files, nat, 
                            require_conversion, conversion_factor):

    f = np.zeros((nat, 3))

    for search_target in xml_files:

        xml = etree.parse(search_target)
        root = xml.getroot()

        for elems in root.findall('calculation/varray'):
            str_force = []
            if elems.get('name') == "forces":
                for elems2 in elems.findall('v'):
                    str_force.append(elems2.text)

            n = len(str_force)

            for i in range(n):
                f[i][:] = [float(t) for t in str_force[i].split()]

            if require_conversion:
                f *= conversion_factor

            for i in range(n):
                print "%15.8E %15.8E %15.8E" % (f[i][0],
                                                f[i][1],
                                                f[i][2])


def print_energies_VASP(xml_files, require_conversion, conversion_factor):

    print "# Etot, Ekin"

    for search_target in xml_files:

        xml = etree.parse(search_target)
        root = xml.getroot()

        for elems in root.findall('calculation/energy'):
            etot = 'N/A'
            ekin = 'N/A'

            for elems2 in elems.findall('i'):
                if elems2.get('name') == "e_fr_energy":
                    etot = elems2.text
                if elems2.get('name') == "kinetic":
                    ekin = elems2.text

            if not require_conversion:
                print "%s %s" % (etot, ekin)
            else:
                if etot != 'N/A':
                    val_etot = float(etot) * conversion_factor
                    print "%15.8E" % val_etot,
                else:
                    print "%s" % etot,

                if ekin != 'N/A':
                    val_ekin = float(ekin) * conversion_factor
                    print "%15.8E" % val_ekin
                else:
                    print "%s" % ekin


# Functions for Quantum-ESPRESSO

def read_original_QE_mod(file_in):

        # Parse general options
    tags = ["ATOMIC_SPECIES", "ATOMIC_POSITIONS", "K_POINTS",
            "CELL_PARAMETERS", "OCCUPATIONS", "CONSTRAINTS", "ATOMIC_FORCES"]

    list_SYSTEM = get_namelist(file_in, "&SYSTEM")
    list_CELL_PARAMETERS = get_options("CELL_PARAMETERS", tags, file_in)
    list_ATOMIC_POSITIONS = get_options("ATOMIC_POSITIONS", tags, file_in)

    ibrav, celldm, nat, ntyp = get_system_info(list_SYSTEM)
    lavec = gen_lattice_vector(ibrav, celldm, list_CELL_PARAMETERS)
    kd_symbol, x0 = get_fractional_coordinate(lavec,
                                              nat,
                                              list_ATOMIC_POSITIONS,
                                              celldm[0])

    return celldm[0], lavec, nat, x0


def print_displacements_QE(pwout_files, alat, lavec, nat, x0,
                           require_conversion, conversion_factor):

    Bohr_to_angstrom = 0.5291772108

    x0 = np.round(x0, 8)
    x = np.zeros((nat, 3))
    disp = np.zeros((nat, 3))

    lavec /= Bohr_to_angstrom
    lavec_transpose = lavec.transpose()
    lavec_transpose_inv = np.linalg.inv(lavec_transpose)

    search_flag = "site n.     atom                  positions (alat units)"
    search_flag2 = "ATOMIC_POSITIONS (crystal)"

    for search_target in pwout_files:

        found_tag = False
        x_list = []
        num_data_disp = 0

        f = open(search_target, 'r')

        line = f.readline()

        while line:

            if search_flag in line:
                found_tag = True

                for i in range(nat):
                    line = f.readline()
                    x[i][:] = [float(t) for t in line.rstrip().split()[6:9]]

                break

            line = f.readline()

        if not found_tag:
            print "%s tag not found in %s" % (search_flag, search_target)
            exit(1)

        x = alat * np.dot(x, lavec_transpose_inv)

        disp = x - x0
        for i in range(nat):
            disp[i][:] = [refold(disp[i][j]) for j in range(3)]

        disp = np.dot(disp, lavec_transpose)

        if require_conversion:
            disp *= conversion_factor

        for i in range(nat):
            print "%15.7F %15.7F %15.7F" % (disp[i][0], 
                                            disp[i][1], 
                                            disp[i][2])

        # Search other entries containing atomis position

        while line:

            if search_flag2 in line:
                num_data_disp += 1

                for i in range(nat):
                    line = f.readline()
                    x[i][:] = [float(t) for t in line.rstrip().split()[1:4]]
                    for j in range(3):
                        x_list.append(x[i][j])

            line = f.readline()
            
        if num_data_disp > 1:
            icount = 0
            for step in range(num_data_disp-1):
                for i in range(nat):
                    for j in range(3):
                        x[i][j] = x_list[icount]
                        icount += 1

                disp = x - x0
                for i in range(nat):
                    disp[i][:] = [refold(disp[i][j]) for j in range(3)]

                disp = np.dot(disp, lavec_transpose)

                if require_conversion:
                    disp *= conversion_factor

                for i in range(nat):
                    print "%15.7F %15.7F %15.7F" % (disp[i][0], 
                                                    disp[i][1], 
                                                    disp[i][2])



def print_atomicforces_QE(str_files, nat, 
                          require_conversion, conversion_factor):

    force = np.zeros((nat, 3))

    search_tag = "Forces acting on atoms (Ry/au):"

    for search_target in str_files:

        found_tag = False

        f = open(search_target, 'r')

        line = f.readline()

        while line:

            if search_tag in line:
                found_tag = True

                f.readline()

                for i in range(nat):
                    line = f.readline()
                    force[i][:] = [float(t) for t in line.rstrip().split()[6:9]]

                if require_conversion:
                    force *= conversion_factor

                for i in range(nat):
                    print "%19.11E %19.11E %19.11E" % (force[i][0],
                                                       force[i][1],
                                                       force[i][2])


            line = f.readline()

        if not found_tag:
            print "%s tag not found in %s" % (search_tag, search_target)
            exit(1)




def print_energies_QE(str_files, require_conversion, conversion_factor):

    search_tag = "!    total energy"

    for search_target in str_files:

        found_tag = False

        with open(search_target) as openfileobject:
            for line in openfileobject:
                if search_tag in line:
                    etot = float(line.rstrip().split()[4])

                    if require_conversion:
                        etot *= conversion_factor

                    print "%19.11E" % etot

                    found_tag = True

        if not found_tag:
            print "%s tag not found in %s" % (search_tag, search_target)
            exit(1)


# Functions for xTAPP

def read_CG_mod(file_in):

    lavec, nat, nkd, list_dummy = read_tappinput(file_in)
    x0, kd, list_dummy = read_atomdata(file_in, nat, nkd)

    return lavec, nat, x0


def print_displacements_xTAPP(str_files, lavec, nat, x0, 
                              require_conversion, conversion_factor):

    Bohr_to_angstrom = 0.5291772108

    x0 = np.round(x0, 8)
    x = np.zeros((nat, 3))

    lavec /= Bohr_to_angstrom
    lavec_transpose = lavec.transpose()

    for search_target in str_files:

        found_tag = False

        f = open(search_target, 'r')

        line = f.readline()

        while line:

            if "atom_position" in line:
                found_tag = True

                for i in range(nat):
                    line = f.readline()
                    x[i, :] = [float(t) for t in line.rstrip().split()[1:]]

                break

            line = f.readline()

        if not found_tag:
            print "atom_position tag not found in %s" % seach_target
            exit(1)

        disp = x - x0
        for i in range(nat):
            disp[i][:] = [refold(disp[i][j]) for j in range(3)]

        disp = np.dot(disp, lavec_transpose)

        if require_conversion:
            disp *= conversion_factor

        for i in range(nat):
            print "%15.7F %15.7F %15.7F" % (disp[i][0], disp[i][1], disp[i][2])


def print_atomicforces_xTAPP(str_files, nat, 
                             require_conversion, conversion_factor):

    force = np.zeros((nat, 3))

    for search_target in str_files:

        found_tag = False

        f = open(search_target, 'r')

        line = f.readline()

        while line:

            if "force" in line:
                found_tag = True

                for i in range(nat):
                    line = f.readline()
                    force[i, :] = [float(t) for t in line.rstrip().split()]

                break

            line = f.readline()

        if not found_tag:
            print "force tag not found in %s" % seach_target
            exit(1)

        if require_conversion:
            force *= conversion_factor

        for i in range(nat):
            print "%19.11E %19.11E %19.11E" % (force[i][0],
                                               force[i][1],
                                               force[i][2])


def print_energies_xTAPP(str_files, require_conversion, conversion_factor):

    for search_target in str_files:

        found_tag = False

        with open(search_target) as openfileobject:
            for line in openfileobject:
                if "total_energy" in line:
                    energy_str = line.rstrip().split()[2]
                    etot = float(energy_str[:-1])

                    if require_conversion:
                        etot *= conversion_factor

                    print "%19.11E" % etot

                    found_tag = True
                    break

        if not found_tag:
            print "total_energy tag not found in %s" % seach_target
            exit(1)


# Other functions

def refold(x):
    if x >= 0.5:
        return x - 1.0
    elif x < -0.5:
        return x + 1.0
    else:
        return x


# Main

if __name__ == "__main__":

    options, args = parser.parse_args()
    file_results = args[0:]

    if len(file_results) == 0:
        print "Usage: extract.py [options] vasprun*.xml \
(or *.pw.out, or *.str)"
        print
        print "For details of available options, please type\n\
$ python displace.py -h"
        exit(1)

    Bohr_radius = 0.52917721092
    Rydberg_to_eV = 13.60569253

    if options.VASP is None and options.QE is None and options.xTAPP is None:
        print "Error : Either --VASP, --QE, or --xTAPP option should be given."
        exit(1)

    elif options.VASP and options.QE or options.VASP and options.xTAPP or\
            options.QE and options.xTAPP:
        print "Error : --VASP, --QE, and --xTAPP \
cannot be given simultaneously."
        exit(1)

    elif options.VASP:
        code = "VASP"
        file_original = options.VASP

        if options.unitname == "eV":
            convert_unit = False
            disp_conv_factor = 1.0
            energy_conv_factor = 1.0

        elif options.unitname == "Rydberg":
            convert_unit = True
            disp_conv_factor = 1.0 / Bohr_radius
            energy_conv_factor = 1.0 / Rydberg_to_eV

        elif options.unitname == "Hartree":
            convert_unit = True
            disp_conv_factor = 1.0 / Bohr_radius
            energy_conv_factor = 0.5 / Rydberg_to_eV

        else:
            print "Invalid options for --unit"
            exit(1)

    elif options.QE:
        code = "QE"
        file_original = options.QE

        if options.unitname == "eV":
            convert_unit = True
            disp_conv_factor = Bohr_radius
            energy_conv_factor = Rydberg_to_eV

        elif options.unitname == "Rydberg":
            convert_unit = False
            disp_conv_factor = 1.0
            energy_conv_factor = 1.0

        elif options.unitname == "Hartree":
            convert_unit = True
            disp_conv_factor = 1.0
            energy_conv_factor = 2.0

        else:
            print "Error : Invalid option for --unit"
            exit(1)

    elif options.xTAPP:
        code = "xTAPP"
        file_original = options.xTAPP

        if options.unitname == "eV":
            convert_unit = True
            disp_conv_factor = Bohr_radius
            energy_conv_factor = 2.0 * Rydberg_to_eV

        elif options.unitname == "Rydberg":
            convert_unit = True
            disp_conv_factor = 1.0
            energy_conv_factor = 2.0

        elif options.unitname == "Hartree":
            convert_unit = False
            disp_conv_factor = 1.0
            energy_conv_factor = 1.0

        else:
            print "Error : Invalid option for --unit"
            exit(1)

    force_conv_factor = energy_conv_factor / disp_conv_factor

    print_disp = False
    print_force = False
    print_energy = False

    if options.get == "disp":
        print_disp = True
    elif options.get == "force":
        print_force = True
    elif options.get == "energy":
        print_energy = True
    else:
        print "Please specify which quantity to extract by the --get option."
        exit(1)

    # Get nat, aa, x_frac0

    if code == "VASP":
        aa, aa_inv, elems, nats, x_frac0 = read_POSCAR(file_original)

    elif code == "QE":
        alat, aa, nat, x_frac0 = read_original_QE_mod(file_original)

    elif code == "xTAPP":
        aa, nat, x_frac0 = read_CG_mod(file_original)

    # Print data

    if print_disp:
        if code == "VASP":
            print_displacements_VASP(file_results, aa, np.sum(nats), x_frac0,
                                     convert_unit, disp_conv_factor)

        elif code == "QE":
            print_displacements_QE(file_results, alat, aa, nat, x_frac0,
                                   convert_unit, disp_conv_factor)

        elif code == "xTAPP":
            print_displacements_xTAPP(file_results, aa, nat, x_frac0,
                                      convert_unit, disp_conv_factor)

    elif print_force:
        if code == "VASP":
            print_atomicforces_VASP(file_results, np.sum(nats),
                                    convert_unit, force_conv_factor)

        elif code == "QE":
            print_atomicforces_QE(file_results, nat,
                                  convert_unit, force_conv_factor)

        elif code == "xTAPP":
            print_atomicforces_xTAPP(file_results, nat,
                                     convert_unit, force_conv_factor)

    elif print_energy:
        if code == "VASP":
            print_energies_VASP(file_results, convert_unit, energy_conv_factor)

        elif code == "QE":
            print_energies_QE(file_results, convert_unit, energy_conv_factor)

        elif code == "xTAPP":
            print_energies_xTAPP(file_results, convert_unit,
                                 energy_conv_factor)
