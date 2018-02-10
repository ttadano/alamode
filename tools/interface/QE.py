#
# QE.py
#
# Interface to Quantum ESPRESSO (http://www.quantum-espresso.org)
#
# Copyright (c) 2014 Terumasa Tadano
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#

import numpy as np


def get_namelist(file_in, namelist_tag):

    list_out = []
    flag_add = False

    with open(file_in) as openfileobject:
        for line in openfileobject:
            line_upper = line.upper()
            if namelist_tag in line_upper:
                flag_add = True
                list_out.append(line)
            elif line.strip() == "/":
                flag_add = False
            elif flag_add:
                list_out.append(line)

    if len(list_out) == 0:
        print("%s field not found" % namelist_tag)
        exit(1)

    list_out.append("/\n")
    return list_out


def gen_lattice_vector(ibrav, celldm, list_CELL_PARAMETERS):
    """.

    Computer lattice vector in units of Angstrom for given ibrav and celldm.
    Doc/INPUT_PW.txt was used as a reference.
    """
    import math

    Bohr_to_angstrom = 0.5291772108

    lavec = np.zeros((3, 3))

    if ibrav == 0:

        if list_CELL_PARAMETERS is None:
            print("CELL_PARAMETERS must be given when ibrav = 0.")
            exit(1)

        else:
            for i in range(3):
                lavec[i][:] = [float(entry) for entry in
                               list_CELL_PARAMETERS[i + 1].rstrip().split()]

            lavec = np.array(lavec)

            mode = list_CELL_PARAMETERS[0].rstrip().split()

            if len(mode) == 1:
                print("Error : Please specify either alat, bohr, or angstrom for CELL_PARAMETERS")
                exit(1)
            else:
                mode_str = mode[1].lower()

            if "alat" in mode_str:

                if not celldm[0]:
                    print("celldm(1) must be given when 'alat' is used for CELL_PARAMETERS")
                    exit(1)

                for i in range(3):
                    for j in range(3):
                        lavec[i][j] *= celldm[0]

            elif "angstrom" in mode_str:

                for i in range(3):
                    for j in range(3):
                        lavec[i][j] /= Bohr_to_angstrom

            elif "bohr" not in mode_str:

                print("Error : Invalid option for CELL_PARAMETERS: %s" % mode[1])
                exit(1)

    elif ibrav == 1:

        if not celldm[0]:
            print("celldm(1) must be given when ibrav = 1.")
            exit(1)

        else:
            a = celldm[0]
            lavec = np.array([[a, 0.0, 0.0],
                              [0.0, a, 0.0],
                              [0.0, 0.0, a]])

    elif ibrav == 2:

        if not celldm[0]:
            print("celldm(1) must be given when ibrav = 2.")
            exit(1)

        else:
            a = celldm[0] / 2.0
            lavec = np.array([[-a, 0.0, a],
                              [0.0, a, a],
                              [-a, a, 0.0]])

    elif ibrav == 3:

        if not celldm[0]:
            print("celldm(1) must be given when ibrav = 3.")
            exit(1)

        else:
            a = celldm[0] / 2.0
            lavec = np.array([[a, a, a],
                              [-a, a, a],
                              [-a, -a, a]])

    elif ibrav == 4:

        if not celldm[0] or not celldm[2]:
            print("celldm(1) and celldm(3) must be given when ibrav = 4.")
            exit(1)

        else:
            a = celldm[0]
            c = celldm[0] * celldm[2]
            lavec = np.array([[a, 0.0, 0.0],
                              [-0.5 * a, math.sqrt(3.) / 2.0 * a, 0.0],
                              [0.0, 0.0, c]])

    elif ibrav == 5 or ibrav == -5:

        if not celldm[0] or not celldm[3]:
            print("celldm(1) and celldm(4) must be given when ibrav = 5, -5.")
            exit(1)

        else:
            a = celldm[0]
            cosalpha = celldm[3]
            tx = a * math.sqrt((1.0 - cosalpha) / 2.)
            ty = a * math.sqrt((1.0 - cosalpha) / 6.)
            tz = a * math.sqrt((1.0 + 2.0 * cosalpha) / 3.)

            if ibrav == 5:
                lavec = np.array([[tx, -ty, tz],
                                  [0.0, 2.0 * ty, tz],
                                  [-tx, -ty, tz]])

            else:
                a_prime = a / math.sqrt(3.0)
                u = tz - 2.0 * math.sqrt(2.0) * ty
                v = tz + math.sqrt(2.0) * ty

                u *= a_prime
                v *= a_prime

                lavec = np.array([[u, v, v],
                                  [v, u, v],
                                  [v, v, u]])

    elif ibrav == 6:

        if not celldm[0] or not celldm[2]:
            print("celldm(1) and celldm(3) must be given when ibrav = 6.")
            exit(1)

        else:
            a = celldm[0]
            c = celldm[0] * celldm[2]
            lavec = np.array([[a, 0.0, 0.0],
                              [0.0, a, 0.0],
                              [0.0, 0.0, c]])

    elif ibrav == 7:

        if not celldm[0] or not celldm[2]:
            print("celldm(1) and celldm(3) must be given when ibrav = 7.")
            exit(1)

        else:
            a = celldm[0]
            c = celldm[0] * celldm[2]
            lavec = np.array([[a / 2.0, -a / 2.0, c / 2.0],
                              [a / 2.0,  a / 2.0, c / 2.0],
                              [-a / 2.0, -a / 2.0, c / 2.0]])

    elif ibrav == 8:

        if not celldm[0] or not celldm[1] or not celldm[2]:
            print("celldm(1), celldm(2), and celldm(3) must be given\
             when ibrav = 8.")
            exit(1)

        else:
            a = celldm[0]
            b = celldm[0] * celldm[1]
            c = celldm[0] * celldm[2]

            lavec = np.array([[a, 0.0, 0.0],
                              [0.0, b, 0.0],
                              [0.0, 0.0, c]])

    elif ibrav == 9 or ibrav == -9:

        if not celldm[0] or not celldm[1] or not celldm[2]:
            print("celldm(1), celldm(2), and celldm(3) must be given\
             when ibrav = 9 or -9.")
            exit(1)

        else:
            a = celldm[0]
            b = celldm[0] * celldm[1]
            c = celldm[0] * celldm[2]

            if ibrav == 9:
                lavec = np.array([[a / 2., b / 2., 0.0],
                                  [-a / 2., b / 2., 0.0],
                                  [0.0, 0.0, c]])
            else:
                lavec = np.array([[a / 2., -b / 2., 0.0],
                                  [a / 2., b / 2., 0.0],
                                  [0.0, 0.0, c]])

    elif ibrav == 10:

        if not celldm[0] or not celldm[1] or not celldm[2]:
            print("celldm(1), celldm(2), and celldm(3) must be given\
             when ibrav = 10.")
            exit(1)

        else:
            a = celldm[0] / 2.0
            b = celldm[0] * celldm[1] / 2.0
            c = celldm[0] * celldm[2] / 2.0
            lavec = np.array([[a, 0.0, c],
                              [a, b, 0.0],
                              [0.0, b, c]])

    elif ibrav == 11:

        if not celldm[0] or not celldm[1] or not celldm[2]:
            print("celldm(1), celldm(2), and celldm(3) must be given\
             when ibrav = 11.")
            exit(1)

        else:
            a = celldm[0] / 2.0
            b = celldm[0] * celldm[1] / 2.0
            c = celldm[0] * celldm[2] / 2.0
            lavec = np.array([[a, b, c],
                              [-a, b, c],
                              [-a, -b, c]])

    elif ibrav == 12:

        if not celldm[0] or not celldm[1] or not celldm[2] or \
           not celldm[3]:
            print("celldm(1), celldm(2), celldm(3), and celldm(4)\
             must be given when ibrav = 12.")
            exit(1)

        else:
            a = celldm[0]
            b = celldm[0] * celldm[1]
            c = celldm[0] * celldm[2]
            gamma = math.acos(celldm[3])
            lavec = np.array([[a, 0.0, 0.0],
                              [b * math.cos(gamma), b * math.sin(gamma), 0.0],
                              [0.0, 0.0, c]])

    elif ibrav == -12:

        if not celldm[0] or not celldm[1] or not celldm[2] or \
           not celldm[4]:
            print("celldm(1), celldm(2), celldm(3), and celldm(5)\
             must be given when ibrav = -12.")
            exit(1)

        else:
            a = celldm[0]
            b = celldm[0] * celldm[1]
            c = celldm[0] * celldm[2]
            beta = math.acos(celldm[4])
            lavec = np.array([[a, 0.0, 0.0],
                              [0.0, b, 0.0],
                              [c * math.cos(beta), 0.0, c * math.sin(beta)]])

    elif ibrav == 13:

        if not celldm[0] or not celldm[1] or not celldm[2] or\
           not celldm[3]:
            print("celldm(1), celldm(2), celldm(3), and celldm(4)\
             must be given when ibrav = 13.")
            exit(1)

        else:
            a = celldm[0]
            b = celldm[0] * celldm[1]
            c = celldm[0] * celldm[2]
            gamma = math.acos(celldm[3])
            lavec = np.array([[a / 2.0, 0.0, -c / 2.0],
                              [b * math.cos(gamma), b * math.sin(gamma), 0.0],
                              [a / 2.0, 0.0, c / 2.0]])

    elif ibrav == 14:

        if not celldm[0] or not celldm[1] or not celldm[2] or \
           not celldm[3] or not celldm[4] or not celldm[5]:
            print("All celldm must be given when ibrav = 14.")
            exit(1)

        else:
            a = celldm[0]
            b = celldm[0] * celldm[1]
            c = celldm[0] * celldm[2]
            alpha = math.acos(celldm[3])
            beta = math.acos(celldm[4])
            gamma = math.acos(celldm[5])

            lavec = np.array([[a, 0.0, 0.0],
                              [b * math.cos(gamma), b * math.sin(gamma), 0.0],
                              [c * math.cos(beta),
                               c * (math.cos(alpha) - math.cos(beta) *
                                    math.cos(gamma)) / math.sin(gamma),
                               c * math.sqrt(1.0 + 2.0 * math.cos(alpha) * math.cos(beta) * math.cos(gamma)
                                             - math.cos(alpha) ** 2 - math.cos(beta) ** 2 - math.cos(gamma) ** 2) / math.sin(gamma)]])

    else:

        print("Invalid ibrav = %s" % ibrav)
        exit(1)

    # Transpose for later use
    lavec = lavec.transpose()

    # Convert to Angstrom unit
    for i in range(3):
        for j in range(3):
            lavec[i][j] *= Bohr_to_angstrom

    return lavec


def get_system_info(list_in):

    list_mod = []

    for obj in list_in:
        obj_split = obj.rstrip().split(',')
        for subobj in obj_split:
            if subobj:
                index = subobj.find('=')
                if index > 0:
                    subobj = subobj[:index] + " = " + subobj[index + 1:]
                list_mod.append(subobj)

    str_input = ""

    for entry in list_mod:
        str_input += entry + " "

    entrylist = str_input.split()

    celldm = [[] for i in range(6)]

    for i in range(len(entrylist)):

        if "ibrav" in entrylist[i]:
            ibrav = int(entrylist[i + 2])

        if "nat" in entrylist[i]:
            nat = int(entrylist[i + 2])

        if "ntyp" in entrylist[i]:
            ntyp = int(entrylist[i + 2])

        if "celldm(1)" in entrylist[i]:
            # Do not assign the value if the comment character '!' 
            # appears in front of the celldm(1) keyword
            has_comment = False
            for elem in list_in:
                if "celldm(1)" in elem:
                    has_comment = ('!' == elem.strip().split()[0][0])

            if not has_comment:
                celldm[0] = float(entrylist[i + 2])

        if "celldm(2)" in entrylist[i]:
            celldm[1] = float(entrylist[i + 2])

        if "celldm(3)" in entrylist[i]:
            celldm[2] = float(entrylist[i + 2])

        if "celldm(4)" in entrylist[i]:
            celldm[3] = float(entrylist[i + 2])

        if "celldm(5)" in entrylist[i]:
            celldm[4] = float(entrylist[i + 2])

        if "celldm(6)" in entrylist[i]:
            celldm[5] = float(entrylist[i + 2])

    return ibrav, celldm, nat, ntyp


def get_options(option_tag, taglists, file_in):

    list_out = []
    flag_add = False

    with open(file_in) as openfileobject:
        for line in openfileobject:

            if option_tag in line:
                flag_add = True
                list_out.append(line)
            elif len(line.split()) > 0 and line.split()[0] in taglists:
                flag_add = False
            elif flag_add:
                list_out.append(line)

    return list_out


def get_fractional_coordinate(aa, N, list_in, a_Bohr):

    Bohr_to_angstrom = 0.5291772108

    list_tmp = list_in[0].rstrip().split()

    if len(list_tmp) == 1:
        print("Error : Please specify either alat, bohr, angstrom, or crystal for ATOMIC_POSITIONS")
        exit(1)
    else:
        mode_str = list_tmp[1].lower()

    if "crystal_sg" in mode_str:
        print("Error : Sorry. 'crystal_sg' is not supported in this script. Please use another option.")
        exit(1)

    xtmp = np.zeros((N, 3))
    kd = []

    for i in range(N):
        list_tmp = list_in[i + 1].rstrip().split()
        kd.append(list_tmp[0])
        xtmp[i][:] = [float(j) for j in list_tmp[1:4]]

    aa_inv = np.linalg.inv(aa)

    if "alat" in mode_str:
        a_angstrom = a_Bohr * Bohr_to_angstrom

        for i in range(3):
            for j in range(3):
                aa_inv[i][j] *= a_angstrom

        for i in range(N):
            xtmp[i][:] = np.dot(xtmp[i][:], aa_inv.transpose())

    elif "bohr" in mode_str:

        for i in range(3):
            for j in range(3):
                aa_inv[i][j] *= Bohr_to_angstrom

        for i in range(N):
            xtmp[i][:] = np.dot(xtmp[i][:], aa_inv.transpose())

    elif "angstrom" in mode_str:

        for i in range(N):
            xtmp[i][:] = np.dot(xtmp[i][:], aa_inv.transpose())

    elif "crystal" not in mode_str:
        print("Error : Invalid option for ATOMIC_POSITIONS: %s" % mode_str)
        exit(1)

    return kd, xtmp


def read_original_QE(file_in):

    # Parse fortran namelists
    list_CONTROL = get_namelist(file_in, "&CONTROL")
    list_SYSTEM = get_namelist(file_in, "&SYSTEM")
    list_ELECTRONS = get_namelist(file_in, "&ELECTRONS")

    # Parse general options
    tags = ["ATOMIC_SPECIES", "ATOMIC_POSITIONS", "K_POINTS",
            "CELL_PARAMETERS", "OCCUPATIONS", "CONSTRAINTS", "ATOMIC_FORCES"]

    list_ATOMIC_SPECIES = get_options("ATOMIC_SPECIES", tags, file_in)
    list_ATOMIC_POSITIONS = get_options("ATOMIC_POSITIONS", tags, file_in)
    list_K_POINTS = get_options("K_POINTS", tags, file_in)
    list_CELL_PARAMETERS = get_options("CELL_PARAMETERS", tags, file_in)
    list_OCCUPATIONS = get_options("OCCUPATIONS", tags, file_in)

    # Get ibrav, celldm, nat, and ntyp
    # and then calculate the lattice vector
    ibrav, celldm, nat, ntyp = get_system_info(list_SYSTEM)
    lavec = gen_lattice_vector(ibrav, celldm, list_CELL_PARAMETERS)
    lavec_inv = np.linalg.inv(lavec)

    # Get fractional coordinate
    kd_symbol, x_frac = get_fractional_coordinate(lavec,
                                                  nat,
                                                  list_ATOMIC_POSITIONS,
                                                  celldm[0])
    list_namelist_merged = []
    list_namelist_merged.extend(list_CONTROL)
    list_namelist_merged.extend(list_SYSTEM)
    list_namelist_merged.extend(list_ELECTRONS)

    return list_namelist_merged, list_ATOMIC_SPECIES, list_K_POINTS, \
        list_CELL_PARAMETERS, list_OCCUPATIONS, \
        nat, lavec, kd_symbol, x_frac, lavec_inv


def generate_QE_input(prefix, suffix, counter, nzerofills, list_namelist,
                      list_ATOMS, list_KP, list_CELL, list_OCCU,
                      nat, kd_symbol, x, u):

    filename = prefix + str(counter).zfill(nzerofills) + "." + suffix
    f = open(filename, 'w')

    for entry in list_namelist:
        f.write(entry)

    for entry in list_ATOMS:
        f.write(entry)

    f.write("ATOMIC_POSITIONS crystal\n")
    for i in range(nat):
        f.write("%s %20.15f %20.15f %20.15f\n" % (kd_symbol[i],
                                                  x[i][0] + u[i, 0],
                                                  x[i][1] + u[i, 1],
                                                  x[i][2] + u[i, 2]))

    for entry in list_KP:
        f.write(entry)
    for entry in list_CELL:
        f.write(entry)
    for entry in list_OCCU:
        f.write(entry)

    f.write("\n")
    f.close()



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


def get_atomicforces_QE(pwout_file, nat):

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
        data0 = get_atomicforces_QE(file_offset, nat)
        try:
            force_offset = np.reshape(data0, (nat, 3))
        except:
            print("File %s contains too many force entries" % file_offset)

    for search_target in str_files:

        force = get_atomicforces_QE(search_target, nat)
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
        print("%s tag not found in %s" % (search_tag, pwout_file))
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


def refold(x):
    if x >= 0.5:
        return x - 1.0
    elif x < -0.5:
        return x + 1.0
    else:
        return x

