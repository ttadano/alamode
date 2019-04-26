#
# OpenMX.py
#
# Interface to OpenMX (http://openmx-square.org)
#
# Copyright (c) 2018 Yuto Tanaka
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#

import numpy as np

"""OpenMX"""
# Function for OpenMX


def read_OpenMX_input(file_original):

    search_target = [
        "atoms.number", "atoms.speciesandcoordinates.unit",
        "<atoms.speciesandcoordinates", "<atoms.unitvectors"
    ]

    # open original file
    f = open(file_original, 'r')

    # set initial patameters
    nat = 0
    lavec_flag = 0
    lavec_row = 0
    lavec = np.zeros([3, 3])

    coord_flag = 0
    coord_row = 0

    # read oroginal file and pull out some infomations
    for line in f:
        ss = line.strip().split()
        # number of atoms
        if len(ss) > 0 and ss[0].lower() == search_target[0]:
            nat = int(ss[1])

        # atomic coordinates
        if coord_flag == 1:
            for j in range(3):
                x_frac0[coord_row][j] = float(ss[j+2])

            coord_row += 1
            if coord_row == nat:
                coord_flag = 0

        # latice vector
        if lavec_flag == 1:
            for i in range(3):
                lavec[lavec_row][i] = float(ss[i])
            lavec_row += 1
            if lavec_row == 3:
                lavec_flag = 0

        # unit of atomic coordinates
        if len(ss) > 0 and ss[0].lower() == search_target[1]:
            coord_unit = ss[1].lower()

        if len(ss) > 0 and ss[0].lower() == search_target[2]:
            coord_flag = 1
            # initialize x_frac0 array
            x_frac0 = np.zeros([nat, 3])

        if len(ss) > 0 and ss[0].lower() == search_target[3]:
            lavec_flag = 1

        if np.linalg.norm(lavec) > 0 and lavec_flag == 0:
            break

    # errors
    if nat == 0:
        print("Could not read dat file properly.")
        exit(1)

    lavec_inv = (np.linalg.inv(lavec)).T
    # convert to frac
    if coord_unit == "ang":
        for i in range(nat):
            x_frac0[i] = np.dot(lavec_inv, x_frac0[i])

    f.close()

    return lavec, lavec_inv, nat, x_frac0


def write_OpenMX_input(prefix, counter, nzerofills, disp, lavec, file_in):

    search_target = [
        "atoms.number", "<atoms.speciesandcoordinates",
        "atoms.speciesandcoordinates.unit", "system.name"
    ]

    filename = prefix + str(counter).zfill(nzerofills) + ".dat"
    fout = open(filename, 'w')
    fin = open(file_in, 'r')

    nat = 0
    coord_flag = 0
    coord_row = 0

    conv = (np.linalg.inv(lavec)).T
    conv_inv = np.linalg.inv(conv)

    for i in range(nat):
        print(np.dot(conv_inv, disp[i]))

    disp[disp < 0] += 1

    for line in fin:
        ss = line.strip().split()
        # number of atoms
        if len(ss) > 0 and ss[0].lower() == search_target[0]:
            nat = int(ss[1])
            x_frac = np.zeros((nat, 3))
            #coord = OrderedDict()
            coord = {}
            for i in range(nat):
                coord[i+1] = []

        # coordinates_unit
        if len(ss) > 0 and ss[0].lower() == search_target[2]:
            coord_unit = ss[1].lower()

        # coordinates
        if coord_flag == 1:
            coord_column = len(ss)
            for i in range(1, coord_column):
                if i > 1:
                    coord[int(ss[0])].append(float(ss[i]))
                else:
                    coord[int(ss[0])].append(ss[i])

            # convert to frac
            if coord_unit == "ang":
                coord[coord_row+1] = np.dot(conv, coord[coord_row+1])

            # add displacement
            for j in range(1, 4):
                coord[coord_row+1][j] += disp[coord_row][j-1]
                coord[coord_row+1][j] = format(coord[coord_row+1][j], '20.16f')

            fout.write(str(coord_row+1) + " ")
            fout.write("  ".join(map(str, coord[coord_row+1])))
            fout.write("\n")
            coord_row += 1
            if coord_row == nat:
                coord_flag = 0

        elif len(ss) > 0 and ss[0].lower() == search_target[3]:
            ss[1] = prefix + str(counter).zfill(nzerofills)
            fout.write("                      ".join(map(str, ss)))
            fout.write("\n")

        else:
            fout.write(line)

        if len(ss) > 0 and ss[0].lower() == search_target[1]:
            coord_flag = 1

    fin.close()
    fout.close()


"""OpenMX"""
# Function for OpenMX


def read_outfile(out_file, nat, column):

    x = np.zeros([nat, 3], dtype=np.float64)

    f = open(out_file, 'r')

    flag = 0
    atom_count = 0
    nat_out = 0
    for line in f:
        ss = line.strip().split()
        if len(ss) > 0:

            if ss[0] == "<coordinates.forces":
                flag = 1
                continue

            if flag == 0:
                continue

            elif flag == 1 and nat_out == 0:
                nat_out = int(ss[0])
                continue

            elif flag == 1 and nat_out > 0:
                for i in range(3):
                    x[atom_count][i] = float(ss[i+column])

                atom_count += 1

        if atom_count == nat:
            break

    f.close()

    return x


# displacements
def get_coordinates_OpenMX(out_file, nat, lavec, conv):
    x = read_outfile(out_file, nat, 2)
    for i in range(nat):
        # convert unit ang to frac
        x[i] = np.dot(conv, x[i])

    return x


def print_displacements_OpenMX(out_files,
                               lavec, lavec_inv, nat, x0,
                               conversion_factor,
                               file_offset):

    vec_refold = np.vectorize(refold)
    lavec_transpose = lavec.transpose()
    conv = lavec_inv
    conv_inv = np.linalg.inv(conv)

    x0 = np.round(x0, 8)

    if file_offset is None:
        disp_offset = np.zeros([nat, 3])
    else:
        x0_offset = get_coordinates_OpenMX(file_offset, nat, lavec, conv)
        try:
            x0_offset = np.reshape(x0_offset, (nat, 3))
        except:
            print("File %s contains too many position entries" % file_offset)
        disp_offset = x0_offset - x0

    for search_target in out_files:

        x = get_coordinates_OpenMX(search_target, nat, lavec, conv)
        #ndata = len(x) / (3 * nat)
        ndata = 1
        #x = np.reshape(x, (1, nat, 3))

        for idata in range(ndata):
            #disp = x[idata, :, :] - x0 - disp_offset
            disp = x - x0 - disp_offset
            disp[disp > 0.96] -= 1.0
            #disp = np.dot(vec_refold(disp), conv_inv)
            for i in range(nat):
                disp[i] = np.dot(conv_inv, disp[i])

            disp[np.absolute(disp) < 1e-5] = 0.0
            disp *= conversion_factor

            for i in range(nat):
                print("%15.7F %15.7F %15.7F" % (disp[i][0],
                                                disp[i][1],
                                                disp[i][2]))


# atomic forces
def get_atomicforces_OpenMX(out_file, nat):
    force = read_outfile(out_file, nat, 5)

    return force


def print_atomicforces_OpenMX(out_files,
                              nat,
                              conversion_factor,
                              file_offset):

    if file_offset is None:
        force_offset = np.zeros((nat, 3))
    else:
        data0 = get_atomicforces_OpenMX(file_offset, nat)
        try:
            force_offset = np.reshape(data0, (nat, 3))
        except:
            print("File %s contains too many force entries" % file_offset)

    for search_target in out_files:
        data = get_atomicforces_OpenMX(search_target, nat)
        #ndata = len(data) / (3 * nat)
        ndata = 1
        #data = np.reshape(data, (ndata, nat, 3))

        for idata in range(ndata):
            #f = data[idata, :, :] - force_offset
            f = data - force_offset

            f *= conversion_factor

            for i in range(nat):
                print("%15.8E %15.8E %15.8E" % (f[i][0],
                                                f[i][1],
                                                f[i][2]))


def print_displacements_and_forces_OpenMX(out_files,
                                          lavec, lavec_inv, nat, x0,
                                          conversion_factor_disp,
                                          conversion_factor_force,
                                          file_offset):

    vec_refold = np.vectorize(refold)
    lavec_transpose = lavec.transpose()
    conv = lavec_inv
    conv_inv = np.linalg.inv(conv)

    x0 = np.round(x0, 8)

    if file_offset is None:
        disp_offset = np.zeros((nat, 3))
        force_offset = np.zeros((nat, 3))
    else:
        x0_offset = get_coordinates_OpenMX(file_offset, nat, lavec, conv)
        force_offset = get_atomicforces_OpenMX(file_offset, nat)

        try:
            x0_offset = np.reshape(x0_offset, (nat, 3))
        except:
            print("File %s contains too many position entries" % file_offset)
        disp_offset = x0_offset - x0

        try:
            force_offset = np.reshape(force_offset, (nat, 3))
        except:
            print("File %s contains too many force entries" % file_offset)

    for search_target in out_files:

        x = get_coordinates_OpenMX(search_target, nat, lavec, conv)
        force = get_atomicforces_OpenMX(search_target, nat)
        ndata = 1

        for idata in range(ndata):
            #disp = x[idata, :, :] - x0 - disp_offset
            disp = x - x0 - disp_offset
            disp[disp > 0.96] -= 1.0
            #disp = np.dot(vec_refold(disp), conv_inv)
            for i in range(nat):
                disp[i] = np.dot(conv_inv, disp[i])

            disp[np.absolute(disp) < 1e-5] = 0.0
            disp *= conversion_factor_disp

            f = force - force_offset
            f *= conversion_factor_force

            for i in range(nat):
                print("%15.7F %15.7F %15.7F %20.8E %15.8E %15.8E" % (disp[i][0],
                                                                     disp[i][1],
                                                                     disp[i][2],
                                                                     f[i][0],
                                                                     f[i][1],
                                                                     f[i][2]))


# total enegy
def get_energies_OpenMX(out_file):

    target = "Utot."
    etot = []

    f = open(out_file, 'r')
    for line in f:
        ss = line.strip().split()
        if len(ss) > 0 and ss[0] == target:
            etot.extend([float(ss[1])])
            break
        else:
            continue

    if len(etot) == 0:
        print("Total energy not found.")
        exit(1)

    return np.array(etot, dtype=np.float)


def print_energies_OpenMX(out_files,
                          conversion_factor,
                          file_offset):

    if file_offset is None:
        etot_offset = 0.0
    else:
        data = get_energies_OpenMX(file_offset)
        if len(data) > 1:
            print("File %s contains too many energy entries" % file_offset)
            exit(1)
        etot_offset = data[0]

    print("# Etot")
    for search_target in out_files:

        etot = get_energies_OpenMX(search_target)

        for idata in range(len(etot)):
            val = etot[idata] - etot_offset

            val *= conversion_factor

            print("%19.11E" % val)


def refold(x):
    if x >= 0.5:
        return x - 1.0
    elif x < -0.5:
        return x + 1.0
    else:
        return x


def get_unit_conversion_factor(str_unit):

    Bohr_radius = 0.52917721067
    Rydberg_to_eV = 13.60569253

    disp_conv_factor = 1.0
    energy_conv_factor = 1.0
    force_conv_factor = 1.0

    if str_unit == "ev":
        disp_conv_factor = 1.0
        energy_conv_factor = 2.0 * Rydberg_to_eV
        force_conv_factor = energy_conv_factor

    elif str_unit == "rydberg":
        disp_conv_factor = 1.0 / Bohr_radius
        energy_conv_factor = 2.0
        force_conv_factor = 2.0

    elif str_unit == "hartree":
        disp_conv_factor = 1.0 / Bohr_radius
        energy_conv_factor = 1.0
        force_conv_factor = 1.0

    else:
        print("This cannot happen")
        exit(1)

    return disp_conv_factor, force_conv_factor, energy_conv_factor


def parse(dat_init, out_files, out_file_offset, str_unit,
          print_disp, print_force, print_energy):

    aa, aa_inv, nat, x_frac0 = read_OpenMX_input(dat_init)

    scale_disp, scale_force, scale_energy = get_unit_conversion_factor(
        str_unit)

    if print_disp is True and print_force is True:
        print_displacements_and_forces_OpenMX(out_files,
                                              aa, aa_inv, nat,
                                              x_frac0,
                                              scale_disp,
                                              scale_force,
                                              out_file_offset)
    elif print_disp is True:
        print_displacements_OpenMX(out_files,
                                   aa, aa_inv, nat,
                                   x_frac0,
                                   scale_disp,
                                   out_file_offset)

    elif print_force is True:
        print_atomicforces_OpenMX(out_files,
                                  nat,
                                  scale_force,
                                  out_file_offset)

    elif print_energy is True:
        print_energies_OpenMX(out_files,
                              scale_energy,
                              out_file_offset)
