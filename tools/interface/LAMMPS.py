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
            if "ITEM:" in line and "ITEM: ATOMS id xu yu zu" not in line:
                add_flag = False
                continue
            elif "ITEM: ATOMS id xu yu zu" in line:
                add_flag = True
                continue

            if add_flag:
                if line.strip():
                    entries = line.strip().split()
                    coord_atom = [int(entries[0]),
                                  [float(t) for t in entries[1:]]]
                    coord.append(coord_atom)

    # This sort is necessary since the order atoms of LAMMPS dump files
    # may change from the input structure file.
    coord_sorted = sorted(coord)
    coord = []
    for coord_atom in coord_sorted:
        coord.extend(coord_atom[1])

    return np.array(coord)


def get_atomicforces_LAMMPS(lammps_dump_file):

    add_flag = False

    force = []

    with open(lammps_dump_file) as f:
        for line in f:
            if "ITEM:" in line and "ITEM: ATOMS id fx fy fz " not in line:
                add_flag = False
                continue
            elif "ITEM: ATOMS id fx fy fz " in line:
                add_flag = True
                continue

            if add_flag:
                if line.strip():
                    entries = line.strip().split()
                    force_atom = [int(entries[0]),
                                  [float(t) for t in entries[1:]]]
                    force.append(force_atom)

    force_sorted = sorted(force)
    force = []
    for force_atom in force_sorted:
        force.extend(force_atom[1])

    return np.array(force)


def get_coordinate_and_force_LAMMPS(lammps_dump_file):

    add_flag = False

    ret = []

    with open(lammps_dump_file) as f:
        for line in f:
            if "ITEM:" in line and "ITEM: ATOMS id xu yu zu fx fy fz" not in line:
                add_flag = False
                continue
            elif "ITEM: ATOMS id xu yu zu fx fy fz" in line:
                add_flag = True
                continue

            if add_flag:
                if line.strip():
                    entries = line.strip().split()
                    data_atom = [int(entries[0]),
                                 [float(t) for t in entries[1:4]],
                                 [float(t) for t in entries[4:]]]

                    ret.append(data_atom)

    # This sort is necessary since the order atoms of LAMMPS dump files
    # may change from the input structure file.
    ret_sorted = sorted(ret)
    ret_x = []
    ret_f = []
    for ret_atom in ret_sorted:
        ret_x.extend(ret_atom[1])
        ret_f.extend(ret_atom[2])

    return np.array(ret_x), np.array(ret_f)


def print_displacements_LAMMPS(lammps_files, nat, x_cart0,
                               conversion_factor, file_offset):

    if file_offset is None:
        disp_offset = np.zeros((nat, 3))
    else:
        _, nat_tmp, x0_offset, _ = read_lammps_structure(file_offset)
        if nat_tmp != nat:
            print("File %s contains too many/few position entries"
                  % file_offset)

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

        # This version supports reading the data from MD trajectory

        for search_target in lammps_files:

            x = get_coordinate_LAMMPS(search_target)
            ndata = len(x) // (3 * nat)
            x = np.reshape(x, (ndata, nat, 3))

            for idata in range(ndata):
                disp = x[idata, :, :] - x_cart0 - disp_offset
                disp *= conversion_factor

                for i in range(nat):
                    print("%20.14f %20.14f %20.14f" % (disp[i, 0],
                                                       disp[i, 1],
                                                       disp[i, 2]))

    else:

        for search_target in lammps_files:

            _, nat_tmp, x_cart, _ = read_lammps_structure(search_target)
            if nat_tmp != nat:
                print("File %s contains too many/few position entries" %
                      search_target)

            disp = x_cart - x_cart0 - disp_offset
            disp *= conversion_factor

            for i in range(nat):
                print("%20.14f %20.14f %20.14f" % (disp[i, 0],
                                                   disp[i, 1],
                                                   disp[i, 2]))


def print_atomicforces_LAMMPS(lammps_files, nat,
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
            f *= conversion_factor

            for i in range(nat):
                print("%19.11E %19.11E %19.11E" % (f[i][0], f[i][1], f[i][2]))


def print_displacements_and_forces_LAMMPS(lammps_files, nat,
                                          x_cart0,
                                          conversion_factor_disp,
                                          conversion_factor_force,
                                          file_offset):

    if file_offset is None:
        disp_offset = np.zeros((nat, 3))
        force_offset = np.zeros((nat, 3))
    else:
        x0_offset, force_offset = get_coordinate_and_force_LAMMPS(file_offset)
        try:
            x0_offset = np.reshape(x0_offset, (nat, 3))
            force_offset = np.reshape(force_offset, (nat, 3))
        except:
            print("File %s contains too many/few entries" % file_offset)

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

        # This version supports reading the data from MD trajectory

        for search_target in lammps_files:

            x, force = get_coordinate_and_force_LAMMPS(search_target)
            ndata = len(x) // (3 * nat)
            x = np.reshape(x, (ndata, nat, 3))
            force = np.reshape(force, (ndata, nat, 3))

            for idata in range(ndata):
                disp = x[idata, :, :] - x_cart0 - disp_offset
                disp *= conversion_factor_disp
                f = force[idata, :, :] - force_offset
                f *= conversion_factor_force

                for i in range(nat):
                    print("%20.14f %20.14f %20.14f %20.8E %15.8E %15.8E" % (disp[i, 0],
                                                                            disp[i, 1],
                                                                            disp[i, 2],
                                                                            f[i, 0],
                                                                            f[i, 1],
                                                                            f[i, 2]))


def get_unit_conversion_factor(str_unit):

    Bohr_radius = 0.52917721067
    Rydberg_to_eV = 13.60569253

    disp_conv_factor = 1.0
    energy_conv_factor = 1.0
    force_conv_factor = 1.0

    if str_unit == "ev":
        disp_conv_factor = 1.0
        energy_conv_factor = 1.0

    elif str_unit == "rydberg":
        disp_conv_factor = 1.0 / Bohr_radius
        energy_conv_factor = 1.0 / Rydberg_to_eV

    elif str_unit == "hartree":
        disp_conv_factor = 1.0 / Bohr_radius
        energy_conv_factor = 0.5 / Rydberg_to_eV

    else:
        print("This cannot happen")
        exit(1)

    force_conv_factor = energy_conv_factor / disp_conv_factor

    return disp_conv_factor, force_conv_factor, energy_conv_factor


def parse(lammps_init, dump_files, dump_file_offset, str_unit,
          print_disp, print_force, print_energy):

    _, nat, x_cart0, _ = read_lammps_structure(lammps_init)
    scale_disp, scale_force, _ = get_unit_conversion_factor(str_unit)

    if print_disp is True and print_force is True:
        print_displacements_and_forces_LAMMPS(dump_files, nat,
                                              x_cart0,
                                              scale_disp,
                                              scale_force,
                                              dump_file_offset)

    elif print_disp is True:
        print_displacements_LAMMPS(dump_files, nat, x_cart0,
                                   scale_disp,
                                   dump_file_offset)

    elif print_force is True:
        print_atomicforces_LAMMPS(dump_files, nat,
                                  scale_force,
                                  dump_file_offset)

    elif print_energy is True:
        print("Error: --get energy is not supported for LAMMPS")
        exit(1)
