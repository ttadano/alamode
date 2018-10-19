#
# xTAPP.py
#
# Interface to xTAPP (http://xtapp.cp.is.s.u-tokyo.ac.jp)
#
# Copyright (c) 2014 Terumasa Tadano
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#

import numpy as np

def read_tappinput(file_in):

    list_tappinput = []
    flag_add = False

    with open(file_in) as openfileobject:
        for line in openfileobject:
            if "main" in line and "data" in line:
                flag_add = True
                list_tappinput.append(line)
            elif "#" in line:
                flag_add = False
            elif flag_add:
                list_tappinput.append(line)

    if len(list_tappinput) == 0:
        print("main data entry not found")
        exit(1)

    list_tappinput_new = []

    for obj in list_tappinput:
        obj_split = obj.rstrip().split(',')
        for subobj in obj_split:
            if subobj:
                list_tappinput_new.append(subobj)

    str_input = ""

    for entry in list_tappinput_new:
        str_input += entry + " "

    entrylist = str_input.split()
    lavec_list = []

    a = 0.0
    nkd = 0
    nat = 0

    # get lattice_factor
    for i in range(len(entrylist)):
        if "lattice_factor" in entrylist[i]:
            a = float(entrylist[i + 2])

        if "lattice_list" in entrylist[i]:
            for j in range(9):
                lavec_list.append(entrylist[i + j + 2])

        if "number_element" in entrylist[i]:
            nkd = int(entrylist[i + 2])

        if "number_atom" in entrylist[i]:
            nat = int(entrylist[i + 2])

    if a == 0.0:
        print("Couldn't read lattice_factor")
        exit(1)
    if nkd == 0:
        print("Couldn't read number_element")
        exit(1)
    if nat == 0:
        print("Couldn't read number_atom")
        exit(1)
    if len(lavec_list) != 9:
        print("Couldn't read lattice_list")
        exit(1)

    lavec = np.zeros((3, 3))

    Bohr_to_angstrom = 0.5291772108
    a *= Bohr_to_angstrom

    for i in range(3):
        for j in range(3):
            lavec[j][i] = a * float(lavec_list[3 * i + j])

    return lavec, nat, nkd, list_tappinput


def read_kpdata(file_in):

    list_kpoint = []
    flag_add = False

    with open(file_in) as openfileobject:
        for line in openfileobject:
            if "k-points" in line.rstrip():
                flag_add = True
                list_kpoint.append(line)
            elif "#" in line.strip():
                flag_add = False
            elif flag_add:
                list_kpoint.append(line)

    if len(list_kpoint) == 0:
        print("k-points data entry not found")
        exit(1)

    return list_kpoint


def read_structure_optimize(file_in):

    list_opt = []
    flag_add = False

    with open(file_in) as openfileobject:
        for line in openfileobject:
            if "struct_opt" in line.rstrip():
                flag_add = True
                list_opt.append(line)
            elif "#" in line.strip():
                flag_add = False
            elif flag_add:
                list_opt.append(line)

    if len(list_opt) == 0:
        print("struct_opt entry not found")
        exit(1)

    list_opt2 = []
    flag_add = False

    with open(file_in) as openfileobject:
        for line in openfileobject:
            if "str_opt_constr" in line.rstrip():
                flag_add = True
                list_opt2.append(line)
            elif "#" in line.strip():
                flag_add = False
            elif flag_add:
                list_opt2.append(line)

    if len(list_opt2) == 0:
        print("str_opt_constr entry not found")
        exit(1)

    return list_opt, list_opt2


def read_atomdata(file_in, nat_in, nkd_in):

    list_atom = []
    flag_add = False

    with open(file_in) as openfileobject:
        for line in openfileobject:
            if "atom" in line and "data" in line:
                flag_add = True
                list_atom.append(line)
            elif "#" in line.strip():
                flag_add = False
            elif flag_add:
                list_atom.append(line)

    if len(list_atom) == 0:
        print("atom data entry not found")
        exit(1)

    x_out = np.zeros((nat_in, 3), dtype=float)
    kd_out = np.zeros(nat_in, dtype=int)

    for i in range(nat_in):
        list_tmp = list_atom[i + nkd_in + 1].rstrip().split()
        kd_out[i] = int(list_tmp[0])
        for j in range(3):
            x_out[i][j] = float(list_tmp[j + 1])

    return x_out, kd_out, list_atom


def read_CG(file_in):

    lavec, nat, nkd, str_tappinput = read_tappinput(file_in)
    str_kpoint = read_kpdata(file_in)
    str_struct_opt, str_opt_constr = read_structure_optimize(file_in)
    x, kd, str_atom = read_atomdata(file_in, nat, nkd)

    str_header = ""

    for entry in str_tappinput:
        str_header += entry
    for entry in str_kpoint:
        str_header += entry
    for entry in str_struct_opt:
        str_header += entry
    for entry in str_opt_constr:
        str_header += entry
    for i in range(nkd + 1):
        str_header += str_atom[i]

    lavec = np.matrix(lavec)
    lavec_inv = np.array(lavec.I)

    return str_header, nat, nkd, lavec, lavec_inv, x, kd


def gen_CG(prefix, suffix, counter, nzerofills, str_header,
           nat, kd, x, u, nsym, symop, denom_tran, has_inv):

    filename = prefix + str(counter).zfill(nzerofills) + "." + suffix
    f = open(filename, 'w')
    f.write("%s" % str_header)

    for i in range(nat):
        f.write("%i %20.15f %20.15f %20.15f\n" % (kd[i], 
                                                  x[i][0] + u[i, 0],
                                                  x[i][1] + u[i, 1],
                                                  x[i][2] + u[i, 2]))

    f.write("# symmetry data\n")
    f.write("&symmetry\n")
    f.write("  number_sym_op = %i\n" % nsym)
    f.write("  has_inversion = %i\n" % has_inv)
    f.write("  denom_trans = %i\n" % denom_tran)
    f.write("/\n")

    mat_tmp = np.zeros((3, 3), dtype=int)

    for elems in symop:
        for i in range(3):
            for j in range(3):
                mat_tmp[i][j] = elems[3 * i + j]

        mat_inv = np.matrix(mat_tmp).I

        for i in range(3):
            for j in range(3):
                f.write("%4i" % mat_inv[i, j])

        f.write("   ")
        for i in range(3):
            f.write("%4i" % elems[9 + i])

        f.write("\n")

    f.write("\n")
    f.close()


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
            disp *= conversion_factor

            for i in range(nat):
                print("%15.7F %15.7F %15.7F" % (disp[i][0],
                                                disp[i][1],
                                                disp[i][2]))


def get_atomicforces_xTAPP(str_file, nat):

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
                             conversion_factor,
                             file_offset):

    if file_offset is None:
        force_offset = np.zeros((nat, 3))
    else:
        data = get_atomicforces_xTAPP(file_offset, nat)
        try:
            force_offset = np.reshape(data, (nat, 3))
        except:
            print("File %s contains too many position entries" % file_offset)

    for search_target in str_files:

        force = get_atomicforces_xTAPP(search_target, nat)
        ndata = len(force) // (3 * nat)
        force = np.reshape(force, (ndata, nat, 3))

        for idata in range(ndata):
            f = force[idata, :, :] - force_offset
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


def print_displacements_and_forces_xTAPP(str_files,
                                         lavec, nat, x0,
                                         conversion_factor_disp,
                                         conversion_factor_force,
                                         file_offset):

    Bohr_to_angstrom = 0.5291772108
    vec_refold = np.vectorize(refold)

    lavec /= Bohr_to_angstrom
    lavec_transpose = lavec.transpose()

    x0 = np.round(x0, 8)

    if file_offset is None:
        disp_offset = np.zeros((nat, 3))
        force_offset = np.zeros((nat, 3))

    else:
        x0_offset = get_coordinates_xTAPP(file_offset, nat)
        force_offset = get_atomicforces_xTAPP(file_offset, nat)
        try:
            x0_offset = np.reshape(x0_offset, (nat, 3))
        except:
            print("File %s contains too many position entries" % file_offset)
        disp_offset = x0_offset - x0

        try:
            force_offset = np.reshape(force_offset, (nat, 3))
        except:
            print("File %s contains too many position entries" % file_offset)

        

    for search_target in str_files:

        x = get_coordinates_xTAPP(search_target, nat)
        force = get_atomicforces_xTAPP(search_target, nat)

        ndata_disp = len(x) // (3 * nat)
        ndata_force = len(force) // (3 * nat)

        if ndata_disp != ndata_force:
            print("Error: The number of entries of displacement and force is inconsistent.")
            print("Ndata disp : %d, Ndata force : %d" % (ndata_disp, ndata_force))
            exit(1)

        ndata = ndata_disp
        x = np.reshape(x, (ndata, nat, 3))
        force = np.reshape(force, (ndata, nat, 3))           

        for idata in range(ndata):
            disp = x[idata, :, :] - x0 - disp_offset
            disp = np.dot(vec_refold(disp), lavec_transpose)
            disp *= conversion_factor
            f = force[idata, :, :] - force_offset
            f *= conversion_factor

            for i in range(nat):
                print("%15.7F %15.7F %15.7F %20.8E %15.8E %15.8E" % (disp[i][0],
                                                                     disp[i][1],
                                                                     disp[i][2],
                                                                     f[i][0],
                                                                     f[i][1],
                                                                     f[i][2]))


def print_energies_xTAPP(str_files,
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
        disp_conv_factor = Bohr_radius
        energy_conv_factor = 2.0 * Rydberg_to_eV

    elif str_unit == "rydberg":
        disp_conv_factor = 1.0
        energy_conv_factor = 2.0

    elif str_unit == "hartree":
        disp_conv_factor = 1.0
        energy_conv_factor = 1.0

    else:
        print("This cannot happen")
        exit(1)
    
    force_conv_factor = energy_conv_factor / disp_conv_factor

    return disp_conv_factor, force_conv_factor, energy_conv_factor

def parse(cg_init, structure_files, structure_file_offset, str_unit,
          print_disp, print_force, print_energy):
    
    aa, nat, x_frac0 = xtapp.read_CG_mod(cg_init)
    
    scale_disp, scale_force, scale_energy = get_unit_conversion_factor(str_unit)

    if print_disp == True and print_force == True:
        print_displacements_and_forces_xTAPP(structure_files,
                                            aa, nat,
                                            x_frac0,
                                            scale_disp,
                                            scale_force,
                                            structure_file_offset)

    elif print_disp == True:
        print_displacements_xTAPP(structure_results, 
                                  aa, nat, x_frac0,
                                  scale_disp, 
                                  structure_file_offset)

    elif print_force == True:
        print_atomicforces_xTAPP(structure_files, 
                                 nat,
                                 scale_force, 
                                 structure_file_offset)

    elif print_energy == True:
        print_energies_xTAPP(structure_files, 
                             scale_energy, 
                             structure_file_offset)

