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

parser.add_option('--get',
                  help="specify which quantity to extract. \
                        Available options are 'disp', 'force' and 'energy'.")


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

# Functions for VASP (https://www.vasp.at)


def get_coordinate_VASP(xml_file, nat):

    x = []

    try:
        xml = etree.parse(xml_file)
        root = xml.getroot()

        for elems in root.findall('calculation/structure/varray'):
            str_coord = [elems2.text for elems2 in elems.findall('v')]
            n = len(str_coord)

            for i in range(n):
                x.extend([t for t in str_coord[i].split()])

        return np.array(x, dtype=np.float)

    except:
        print("Error in reading atomic positions from the XML file: %s" % xml_file)


def print_displacements_VASP(xml_files,
                             lavec, nat, x0,
                             require_conversion,
                             conversion_factor,
                             file_offset):

    x0 = np.round(x0, 8)
    lavec_transpose = lavec.transpose()
    vec_refold = np.vectorize(refold)

    if file_offset is None:
        disp_offset = np.zeros((nat, 3))
    else:
        x0_offset = get_coordinate_VASP(file_offset, nat)
        try:
            x0_offset = np.reshape(x0_offset, (nat, 3))
        except:
            print("File %s contains too many position entries" % file_offset)
        disp_offset = x0_offset - x0

    for search_target in xml_files:

        x = get_coordinate_VASP(search_target, nat)
        ndata = len(x) // (3 * nat)
        x = np.reshape(x, (ndata, nat, 3))

        for idata in range(ndata):
            disp = x[idata, :, :] - x0 - disp_offset
            disp = np.dot(vec_refold(disp), lavec_transpose)

            if require_conversion:
                disp *= conversion_factor

            for i in range(nat):
                print("%15.7F %15.7F %15.7F" % (disp[i, 0],
                                                disp[i, 1],
                                                disp[i, 2]))


def get_atomicforces_VASP(xml_file):

    f = []

    try:
        xml = etree.parse(xml_file)
        root = xml.getroot()

        for elems in root.findall('calculation/varray'):
            if elems.get('name') == "forces":
                str_force = [elems2.text for elems2 in elems.findall('v')]

                for i in range(len(str_force)):
                    f.extend([t for t in str_force[i].split()])
                
        return np.array(f, dtype=np.float)

    except:
        print("Error in reading atomic forces from the XML file: %s" % xml_file)


def print_atomicforces_VASP(xml_files,
                            nat,
                            require_conversion,
                            conversion_factor,
                            file_offset):

    if file_offset is None:
        force_offset = np.zeros((nat, 3))
    else:
        data0 = get_atomicforces_VASP(file_offset)
        try:
            force_offset = np.reshape(data0, (nat, 3))
        except:
            print("File %s contains too many force entries" % file_offset)

    for search_target in xml_files:

        data = get_atomicforces_VASP(search_target)
        ndata = len(data) // (3 * nat)
        data = np.reshape(data, (ndata, nat, 3))

        for idata in range(ndata):
            f = data[idata, :, :] - force_offset

            if require_conversion:
                f *= conversion_factor

            for i in range(nat):
                print("%15.8E %15.8E %15.8E" % (f[i][0],
                                                f[i][1],
                                                f[i][2]))


def get_energies_VASP(xml_file):

    etot_array = []
    ekin_array = []

    try:
        xml = etree.parse(xml_file)
        root = xml.getroot()

        for elems in root.findall('calculation/energy'):
            etot = 'N/A'
            ekin = 'N/A'

            for elems2 in elems.findall('i'):
                if elems2.get('name') == "e_fr_energy":
                    etot = elems2.text
                if elems2.get('name') == "kinetic":
                    ekin = elems2.text

            etot_array.append(etot)
            ekin_array.append(ekin)

        return etot_array, ekin_array
    except:
        print("Error in reading energies from the XML file: %s" % xml_file)


def print_energies_VASP(xml_files,
                        require_conversion,
                        conversion_factor,
                        file_offset):

    print("# Etot, Ekin")

    etot_offset = 0.0
    ekin_offset = 0.0

    if file_offset:
        etot, ekin = get_energies_VASP(file_offset)
        if len(etot) > 1 or len(ekin) > 1:
            print("File %s contains too many energy entries" % file_offset)
            exit(1)
        if etot[0] != 'N/A':
            etot_offset = float(etot[0])
        if ekin[0] != 'N/A':
            ekin_offset = float(ekin[0])

    for search_target in xml_files:

        etot, ekin = get_energies_VASP(search_target)

        for i in range(len(etot)):
            if etot[i] != 'N/A':
                val_etot = float(etot[i]) - etot_offset
                if require_conversion:
                    print("%15.8E" % (val_etot * conversion_factor), end=' ')
                else:
                    print("%15.8E" % val_etot, end=' ')
            else:
                print("%s" % etot[i], end=' ')

            if ekin[i] != 'N/A':
                val_ekin = float(ekin[i]) - ekin_offset
                if require_conversion:
                    print("%15.8E" % (val_ekin * conversion_factor))
                else:
                    print("%15.8E" % val_ekin)
            else:
                print("%s" % ekin[i])

# end functions for VASP

# Functions for Quantum-ESPRESSO (http://www.quantum-espresso.org)


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


def get_coordinates_QE(pwout_file, nat):

    search_flag = "site n.     atom                  positions (alat units)"
    search_flag2 = "ATOMIC_POSITIONS (crystal)"

    x = np.zeros((nat, 3))

    num_data_disp = 0
    basis = ""
    found_tag = False

    f = open(pwout_file, 'r')
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
        print("%s tag not found in %s" % (search_flag, pwout_file))
        exit(1)

    x_additional = []

    # Search other entries containing atomic position
    while line:

        if search_flag2 in line:

            if not basis:
                basis = line.rstrip().split()[1]

            num_data_disp += 1

            for i in range(nat):
                line = f.readline()
                x_additional.extend([t for t in line.rstrip().split()[1:4]])

        line = f.readline()

    f.close()

    return x, np.array(x_additional, dtype=np.float), num_data_disp, basis


def print_displacements_QE(pwout_files,
                           alat, lavec, nat, x0,
                           require_conversion,
                           conversion_factor,
                           file_offset):

    import math
    Bohr_to_angstrom = 0.5291772108
    vec_refold = np.vectorize(refold)

    x0 = np.round(x0, 8)

    lavec /= Bohr_to_angstrom
    lavec_transpose = lavec.transpose()
    lavec_transpose_inv = np.linalg.inv(lavec_transpose)

    if not alat:
        # if celldm[0] is empty, calculate it from lattice vector
        alat = math.sqrt(np.dot(lavec_transpose[0][:], lavec_transpose[0][:]))

    if file_offset is None:
        disp_offset = np.zeros((nat, 3))
    else:
        x_offset, x_tmp, ndata_offset, basis_tmp = get_coordinates_QE(
            file_offset, nat)
        if ndata_offset > 1:
            print("File %s contains too many position entries" % file_offset)
            exit(1)
        else:
            x_offset = alat * np.dot(x_offset, lavec_transpose_inv)
            disp_offset = x_offset - x0

    for search_target in pwout_files:

        x, x_additional, num_data_disp, basis = get_coordinates_QE(
            search_target, nat)
        x = alat * np.dot(x, lavec_transpose_inv)

        disp = x - x0 - disp_offset
        disp = np.dot(vec_refold(disp), lavec_transpose)

        if require_conversion:
            disp *= conversion_factor

        for i in range(nat):
            print("%15.7F %15.7F %15.7F" % (disp[i][0],
                                            disp[i][1],
                                            disp[i][2]))

        if num_data_disp > 1:

            if "alat" in basis:
                conversion_mat = alat * lavec_transpose_inv
            elif "bohr" in basis:
                conversion_mat = lavec_transpose_inv
            elif "angstrom" in basis:
                conversion_mat = lavec_transpose_inv / Bohr_to_angstrom
            elif "crystal" in basis:
                conversion_mat = np.identity(3)
            else:
                print("This cannot happen.")
                exit(1)

            x_additional = np.reshape(x_additional, (num_data_disp, nat, 3))

            for step in range(num_data_disp - 1):
                x = x_additional[step, :, :]
                x = np.dot(x, conversion_mat)
                disp = x - x0 - disp_offset
                disp = np.dot(vec_refold(disp), lavec_transpose)

                if require_conversion:
                    disp *= conversion_factor

                for i in range(nat):
                    print("%15.7F %15.7F %15.7F" % (disp[i][0],
                                                    disp[i][1],
                                                    disp[i][2]))


def get_atomicforces_QE(pwout_file):

    search_tag = "Forces acting on atoms (Ry/au):"
    search_tag_QE6 = "Forces acting on atoms (cartesian axes, Ry/au):"

    found_tag = False

    f = open(pwout_file, 'r')
    line = f.readline()

    force = []

    while line:

        if search_tag in line or search_tag_QE6 in line:
            found_tag = True

            f.readline()

            for i in range(nat):
                line = f.readline()
                force.extend([t for t in line.rstrip().split()[6:9]])

        line = f.readline()

    f.close()

    if not found_tag:
        print("following search tags not found in %s" %  pwout_file)
        print(search_tag)
        print(search_tag_QE6)
        exit(1)

    return np.array(force, dtype=np.float)


def print_atomicforces_QE(str_files,
                          nat,
                          require_conversion,
                          conversion_factor,
                          file_offset):

    if file_offset is None:
        force_offset = np.zeros((nat, 3))
    else:
        data0 = get_atomicforces_QE(file_offset)
        try:
            force_offset = np.reshape(data0, (nat, 3))
        except:
            print("File %s contains too many force entries" % file_offset)

    for search_target in str_files:

        force = get_atomicforces_QE(search_target)
        ndata = len(force) // (3 * nat)
        force = np.reshape(force, (ndata, nat, 3))

        for idata in range(ndata):
            f = force[idata, :, :] - force_offset

            if require_conversion:
                f *= conversion_factor

            for i in range(nat):
                print("%19.11E %19.11E %19.11E" % (f[i][0],
                                                   f[i][1],
                                                   f[i][2]))


def get_energies_QE(pwout_file):

    search_tag = "!    total energy"

    found_tag = False

    etot = []

    with open(pwout_file) as openfileobject:
        for line in openfileobject:
            if search_tag in line:
                etot.extend([line.rstrip().split()[4]])
                found_tag = True

    if not found_tag:
        print("%s tag not found in %s" % (search_tag, search_target))
        exit(1)

    return np.array(etot, dtype=np.float)


def print_energies_QE(str_files,
                      require_conversion,
                      conversion_factor,
                      file_offset):

    if file_offset is None:
        etot_offset = 0.0
    else:
        data = get_energies_QE(file_offset)
        if len(data) > 1:
            print("File %s contains too many energy entries" % file_offset)
            exit(1)
        etot_offset = data[0]

    print("# Etot")
    for search_target in str_files:

        etot = get_energies_QE(search_target)

        for idata in range(len(etot)):
            val = etot[idata] - etot_offset

            if require_conversion:
                val *= conversion_factor

            print("%19.11E" % val)

# end functions for QE

# Functions for xTAPP (http://xtapp.cp.is.s.u-tokyo.ac.jp)


def read_CG_mod(file_in):

    lavec, nat, nkd, list_dummy = read_tappinput(file_in)
    x0, kd, list_dummy = read_atomdata(file_in, nat, nkd)

    return lavec, nat, x0


def get_coordinates_xTAPP(str_file, nat):

    found_tag = False
    f = open(str_file, 'r')
    line = f.readline()

    x = []

    while line:

        if "atom_position" in line:
            found_tag = True

            for i in range(nat):
                line = f.readline()
                x.extend([t for t in line.rstrip().split()[1:]])

            break

        line = f.readline()

    if not found_tag:
        print("atom_position tag not found in %s" % str_file)
        exit(1)

    f.close()

    return np.array(x, dtype=np.float)


def print_displacements_xTAPP(str_files,
                              lavec, nat, x0,
                              require_conversion,
                              conversion_factor,
                              file_offset):

    Bohr_to_angstrom = 0.5291772108
    vec_refold = np.vectorize(refold)

    lavec /= Bohr_to_angstrom
    lavec_transpose = lavec.transpose()

    x0 = np.round(x0, 8)

    if file_offset is None:
        disp_offset = np.zeros((nat, 3))
    else:
        x0_offset = get_coordinates_xTAPP(file_offset, nat)
        try:
            x0_offset = np.reshape(x0_offset, (nat, 3))
        except:
            print("File %s contains too many position entries" % file_offset)
        disp_offset = x0_offset - x0

    for search_target in str_files:

        x = get_coordinates_xTAPP(search_target, nat)
        ndata = len(x) // (3 * nat)
        x = np.reshape(x, (ndata, nat, 3))

        for idata in range(ndata):
            disp = x[idata, :, :] - x0 - disp_offset
            disp = np.dot(vec_refold(disp), lavec_transpose)

            if require_conversion:
                disp *= conversion_factor

            for i in range(nat):
                print("%15.7F %15.7F %15.7F" % (disp[i][0],
                                                disp[i][1],
                                                disp[i][2]))


def get_atomicforces_xTAPP(str_file):

    found_tag = False

    f = open(str_file, 'r')
    line = f.readline()

    force = []

    while line:

        if "force" in line:
            found_tag = True

            for i in range(nat):
                line = f.readline()
                force.extend([t for t in line.rstrip().split()])

            break

        line = f.readline()

    if not found_tag:
        print("force tag not found in %s" % str_file)
        exit(1)

    f.close()

    return np.array(force, dtype=np.float)


def print_atomicforces_xTAPP(str_files,
                             nat,
                             require_conversion,
                             conversion_factor,
                             file_offset):

    if file_offset is None:
        force_offset = np.zeros((nat, 3))
    else:
        data = get_atomicforces_xTAPP(file_offset)
        try:
            force_offset = np.reshape(data, (nat, 3))
        except:
            print("File %s contains too many position entries" % file_offset)

    for search_target in str_files:

        force = get_atomicforces_xTAPP(search_target)
        ndata = len(force) // (3 * nat)
        force = np.reshape(force, (ndata, nat, 3))

        for idata in range(ndata):
            f = force[idata, :, :] - force_offset

            if require_conversion:
                f *= conversion_factor

                for i in range(nat):
                    print("%19.11E %19.11E %19.11E" % (f[i][0],
                                                       f[i][1],
                                                       f[i][2]))


def get_energies_xTAPP(str_file):

    search_tag = "total_energy"

    found_tag = False

    etot = []

    with open(str_file) as openfileobject:
        for line in openfileobject:
            if search_tag in line:
                energy_str = line.rstrip().split()[2]
                etot.extend([energy_str[:-1]])
                found_tag = True

    if not found_tag:
        print("%s tag not found in %s" % (search_tag, str_file))
        exit(1)

    return np.array(etot, dtype=np.float)


def print_energies_xTAPP(str_files,
                         require_conversion,
                         conversion_factor,
                         file_offset):

    if file_offset is None:
        etot_offset = 0.0
    else:
        data = get_energies_xTAPP(file_offset)
        if len(data) > 1:
            print("File %s contains too many energy entries" % file_offset)
            exit(1)
        etot_offset = data[0]

    print("# Etot")
    for search_target in str_files:

        etot = get_energies_xTAPP(search_target)

        for idata in range(len(etot)):
            val = etot[idata] - etot_offset

            if require_conversion:
                val *= conversion_factor

            print("%19.11E" % val)


# end functions for xTAPP


# Functions for LAMMPS


def get_coordinate_LAMMPS(lammps_dump_file):
    
    add_flag = False

    coord = []

    with open(lammps_dump_file) as f:
        for line in f:
            if "ITEM:" in line and "ITEM: ATOMS xu yu zu" not in line:
                add_flag = False
                continue
            elif "ITEM: ATOMS xu yu zu" in line:
                add_flag = True
                continue

            if add_flag:
                if line.strip():
                    coord.extend([float(t) for t in line.strip().split()])

    return np.array(coord)


def get_atomicforces_LAMMPS(lammps_dump_file):
    
    add_flag = False

    force = []

    with open(lammps_dump_file) as f:
        for line in f:
            if "ITEM:" in line and "ITEM: ATOMS fx fy fz " not in line:
                add_flag = False
                continue
            elif "ITEM: ATOMS fx fy fz " in line:
                add_flag = True
                continue

            if add_flag:
                if line.strip():
                    force.extend([float(t) for t in line.strip().split()])

    return np.array(force)


def print_displacements_LAMMPS(lammps_files,
                               nat, x_cart0,
                               require_conversion,
                               conversion_factor,
                               file_offset):

    if file_offset is None:
        disp_offset = np.zeros((nat, 3))
    else:
        dummy, nat_tmp, x0_offset, kd_offset = read_lammps_structure(file_offset)
        if nat_tmp != nat:
            print("File %s contains too many/few position entries" % file_offset)

        disp_offset = x0_offset - x_cart0

    # Automatic detection of the input format 

    is_dumped_file = False
    f = open(lammps_files[0], 'r')
    for line in f:
        if "ITEM: TIMESTEP" in line:
            is_dumped_file = True
            break
    f.close()

    if is_dumped_file:
        
        ## This version supports reading the data from MD trajectory
        
         for search_target in lammps_files:
        
            x = get_coordinate_LAMMPS(search_target)
            ndata = len(x) // (3 * nat)
            x = np.reshape(x, (ndata, nat, 3))

            for idata in range(ndata):
                disp = x[idata, :, :] - x_cart0 - disp_offset

                if require_conversion:
                    disp *= conversion_factor

                for i in range(nat):
                    print("%20.14f %20.14f %20.14f" % (disp[i, 0],
                                                       disp[i, 1],
                                                       disp[i, 2]))
        
    else:

        for search_target in lammps_files:

            dummy, nat_tmp, x_cart, kd_tmp = read_lammps_structure(search_target)
            if nat_tmp != nat:
                print("File %s contains too many/few position entries" % search_target)

            disp = x_cart - x_cart0 - disp_offset

            if require_conversion:
                disp *= conversion_factor

            for i in range(nat):
                print("%20.14f %20.14f %20.14f" % (disp[i, 0],
                                                   disp[i, 1],
                                                   disp[i, 2]))


def print_atomicforces_LAMMPS(lammps_files, nat, 
                              require_conversion,
                              conversion_factor,
                              file_offset):

    if file_offset is None:
        force_offset = np.zeros((nat, 3))
    else:
        data = get_atomicforces_LAMMPS(file_offset)
        try:
            force_offset = np.reshape(data, (nat, 3))
        except:
            print("File %s contains too many position entries" % file_offset)


    # Automatic detection of the input format 

    is_dumped_file = False
    f = open(lammps_files[0], 'r')
    for line in f:
        if "ITEM: TIMESTEP" in line:
            is_dumped_file = True
            break
    f.close()

    for search_target in lammps_files:
    
        force = get_atomicforces_LAMMPS(search_target)
        ndata = len(force) // (3 * nat)
        force = np.reshape(force, (ndata, nat, 3))

        for idata in range(ndata):
            f = force[idata, :, :] - force_offset

            if require_conversion:
                f *= conversion_factor

            for i in range(nat):
                print("%19.11E %19.11E %19.11E" % (f[i][0],
                                                    f[i][1],
                                                    f[i][2]))


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
        print("Usage: extract.py [options] vasprun*.xml \
(or *.pw.out, or *.str)")
        print()
        print("For details of available options, please type\n\
$ python displace.py -h")
        exit(1)

    Bohr_radius = 0.52917721067
    Rydberg_to_eV = 13.60569253

    conditions = [options.VASP is None,
                  options.QE is None,
                  options.xTAPP is None,
                  options.LAMMPS is None]

    if conditions.count(True) == len(conditions):
        print("Error : Either --VASP, --QE, --xTAPP, --LAMMPS option must be given.")
        exit(1)

    elif len(conditions) - conditions.count(True) > 1:
        print("Error : --VASP, --QE, --xTAPP, and --LAMMPS cannot be given simultaneously.")
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
            print("Invalid options for --unit")
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
            energy_conv_factor = 0.5

        else:
            print("Error : Invalid option for --unit")
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
            print("Error : Invalid option for --unit")
            exit(1)

    elif options.LAMMPS:
        code = "LAMMPS"
        file_original = options.LAMMPS

        if options.get == "energy":
            print("--get energy is not supported for LAMMPS")
            exit(1)

        if options.unitname == "eV":
            convert_unit = False
            disp_conv_factor = 1.0
            energy_conv_factor = 1.0

        else:
            convert_unit = True
            disp_conv_factor = 1.0 / Bohr_radius
            energy_conv_factor = 1.0 / Rydberg_to_eV

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
        print("Please specify which quantity to extract by the --get option.")
        exit(1)

    # Get nat, aa, x_frac0

    if code == "VASP":
        aa, aa_inv, elems, nats, x_frac0 = read_POSCAR(file_original)

    elif code == "QE":
        alat, aa, nat, x_frac0 = read_original_QE_mod(file_original)

    elif code == "xTAPP":
        aa, nat, x_frac0 = read_CG_mod(file_original)
    
    elif code == "LAMMPS":
        common_settings, nat, x_cart0, kd = read_lammps_structure(file_original)

    # Print data

    if print_disp:
        if code == "VASP":
            print_displacements_VASP(file_results, aa, np.sum(nats), x_frac0,
                                     convert_unit, disp_conv_factor, options.offset)

        elif code == "QE":
            print_displacements_QE(file_results, alat, aa, nat, x_frac0,
                                   convert_unit, disp_conv_factor, options.offset)

        elif code == "xTAPP":
            print_displacements_xTAPP(file_results, aa, nat, x_frac0,
                                      convert_unit, disp_conv_factor, options.offset)

        elif code == "LAMMPS":
            print_displacements_LAMMPS(file_results, nat, x_cart0,
                                       convert_unit, disp_conv_factor, options.offset)

    elif print_force:
        if code == "VASP":
            print_atomicforces_VASP(file_results, np.sum(nats),
                                    convert_unit, force_conv_factor, options.offset)

        elif code == "QE":
            print_atomicforces_QE(file_results, nat,
                                  convert_unit, force_conv_factor, options.offset)

        elif code == "xTAPP":
            print_atomicforces_xTAPP(file_results, nat,
                                     convert_unit, force_conv_factor, options.offset)

        elif code == "LAMMPS":
            print_atomicforces_LAMMPS(file_results, nat, 
                                      convert_unit, force_conv_factor, options.offset)

    elif print_energy:
        if code == "VASP":
            print_energies_VASP(file_results, convert_unit,
                                energy_conv_factor, options.offset)

        elif code == "QE":
            print_energies_QE(file_results, convert_unit,
                              energy_conv_factor, options.offset)

        elif code == "xTAPP":
            print_energies_xTAPP(file_results, convert_unit,
                                 energy_conv_factor, options.offset)
