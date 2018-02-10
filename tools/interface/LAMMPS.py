#
# LAMMPS.py
#
# Interface to LAMMPS (http://lammps.sandia.gov)
#
# Copyright (c) 2017 Terumasa Tadano
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#

import numpy as np


def read_lammps_structure(file_in):

    f = open(file_in, 'r')
    header_comment = f.readline()

    common_settings = []

    for line in f:
        if "Atoms" in line:
            break
        common_settings.append(line.rstrip())

    atoms = []
    for line in f:
        if line.strip():
            atoms.append(line.rstrip().split())

    atoms = np.array(atoms)
    nat = len(atoms)
    kd = np.array(atoms[:, 1], dtype=np.int)
    x = np.array(atoms[:, 2:5], dtype=np.float64)

    return common_settings, nat, x, kd


def write_lammps_structure(prefix, counter, header, nzerofills,
                           common_settings, nat, kd, x_cart, disp):

    filename = prefix + str(counter).zfill(nzerofills) + ".lammps"
    f = open(filename, 'w')
    f.write("%s\n" % header)

    for line in common_settings:
        f.write("%s\n" % line)

    f.write("%s\n\n" % "Atoms")
    for i in range(nat):
        f.write("%5d %3d" % (i + 1, kd[i]))
        for j in range(3):
            f.write("%20.15f" % (x_cart[i][j] + disp[i][j]))
        f.write("\n")
    f.write("\n")
    f.close()


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


def print_displacements_LAMMPS(lammps_files, nat, x_cart0, require_conversion,
                               conversion_factor, file_offset):

    if file_offset is None:
        disp_offset = np.zeros((nat, 3))
    else:
        dummy, nat_tmp, x0_offset, kd_offset = read_lammps_structure(
            file_offset)
        if nat_tmp != nat:
            print(
                "File %s contains too many/few position entries" % file_offset)

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

            dummy, nat_tmp, x_cart, kd_tmp = read_lammps_structure(
                search_target)
            if nat_tmp != nat:
                print("File %s contains too many/few position entries" %
                      search_target)

            disp = x_cart - x_cart0 - disp_offset

            if require_conversion:
                disp *= conversion_factor

            for i in range(nat):
                print("%20.14f %20.14f %20.14f" % (disp[i, 0], 
                                                   disp[i, 1],
                                                   disp[i, 2]))


def print_atomicforces_LAMMPS(lammps_files, nat, require_conversion,
                              conversion_factor, file_offset):

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
                print("%19.11E %19.11E %19.11E" % (f[i][0], f[i][1], f[i][2]))
