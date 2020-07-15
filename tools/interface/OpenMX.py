#
# OpenMX.py
#
# Interface to OpenMX (http://openmx-square.org)
#
# Copyright (c) 2018 Yuto Tanaka
# Copyright (c) 2020 Terumasa Tadano
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#

import numpy as np


class OpenmxParser(object):

    def __init__(self):
        self._prefix = None
        self._lattice_vector = None
        self._inverse_lattice_vector = None
        self._kd = None
        self._initial_charges = None
        self._common_settings = None
        self._nat = 0
        self._x_fractional = None
        self._counter = 1
        self._nzerofills = 0
        self._disp_conversion_factor = 1.0
        self._energy_conversion_factor = 1.0
        self._force_conversion_factor = 1.0
        self._initial_structure_loaded = False
        self._print_disp = True
        self._print_force = True
        self._print_energy = False
        self._print_born = False
        self._BOHR_TO_ANGSTROM = 0.5291772108
        self._RYDBERG_TO_EV = 13.60569253

    def load_initial_structure(self, file_original):

        search_target = [
            "atoms.number", "atoms.speciesandcoordinates.unit",
            "<atoms.speciesandcoordinates", "<atoms.unitvectors",
            "atoms.speciesandcoordinates>"
        ]

        nat = None
        lavec = []
        common_settings = []
        kd = []
        x_frac0 = []
        initial_charges = []

        # read original file and pull out some information

        with open(file_original, 'r') as f:
            lines = f.read().splitlines()

            for i, line in enumerate(lines):

                if search_target[0] in line.lower():
                    nat = int(line.strip().split()[-1])

                elif search_target[1] in line.lower():
                    coord_unit = line.strip().split()[-1].lower()

                elif search_target[2] in line.lower():
                    ipos_coord = i + 1

                elif search_target[4] in line.lower():
                    fpos_coord = i

                elif search_target[3] in line.lower():
                    ipos_lavec = i + 1

            if nat is None:
                raise RuntimeError("Failed to extract the Atoms.Number value from the file.")

            if nat != (fpos_coord - ipos_coord):
                raise RuntimeError("The number of entries in Atoms.SpeciesAndCoordinates does not match"
                                   "with the Atoms.Number value.")

            for line in lines[ipos_coord:fpos_coord]:
                line_split = line.strip().split()
                kd.append(line_split[1])
                x_frac0.append([float(t) for t in line_split[2:5]])
                initial_charges.append([float(t) for t in line_split[5:7]])

            for line in lines[ipos_lavec:ipos_lavec+3]:
                lavec.append([float(t) for t in line.strip().split()])

            common_settings.append(lines[:ipos_coord])
            common_settings.append(line[fpos_coord:])

        x_frac0 = np.array(x_frac0)
        lavec = np.array(lavec).transpose()
        lavec_inv = np.linalg.inv(lavec)
        initial_charges = np.array(initial_charges)

        # convert to frac
        if coord_unit == "ang":
            for i in range(nat):
                x_frac0[i] = np.dot(x_frac0[i], lavec_inv)

        self._lattice_vector = lavec
        self._inverse_lattice_vector = np.linalg.inv(lavec)
        self._nat = nat
        self._x_fractional = x_frac0
        self._kd = kd
        self._initial_charges = initial_charges
        self._common_settings = common_settings
        self._initial_structure_loaded = True


    def generate_structures(self, prefix, header_list, disp_list):

        self._set_number_of_zerofill(len(disp_list))
        self._prefix = prefix

        for header, disp in zip(header_list, disp_list):
            self._generate_input(header, disp)

    def parse(self, initial_pwin, pwout_files, pwout_file_offset, str_unit,
              output_flags, filter_emin=None, filter_emax=None):

        if not self._initial_structure_loaded:
            self.load_initial_structure(initial_pwin)

        self._set_unit_conversion_factor(str_unit)
        self._set_output_flags(output_flags)

        if self._print_disp or self._print_force:
            self._print_displacements_and_forces(pwout_files,
                                                 pwout_file_offset,
                                                 filter_emin,
                                                 filter_emax)
        elif self._print_energy:
            self._print_energies(pwout_files, pwout_file_offset)




def generate_input(prefix, counter, disp, nzerofills, params_orig, file_in):

    search_target = [
        "atoms.number", "<atoms.speciesandcoordinates",
        "atoms.speciesandcoordinates.unit", "system.name"
    ]

    lavec = params_orig['lavec']

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
            raise RuntimeError("File %s contains too many position entries" % file_offset)
        disp_offset = x0_offset - x0

    for search_target in out_files:

        x = get_coordinates_OpenMX(search_target, nat, lavec, conv)
        ndata = 1

        for idata in range(ndata):
            disp = x - x0 - disp_offset
            disp[disp > 0.96] -= 1.0
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
            raise RuntimeError("File %s contains too many force entries" % file_offset)

    for search_target in out_files:
        data = get_atomicforces_OpenMX(search_target, nat)
        ndata = 1

        for idata in range(ndata):
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
            disp = x - x0 - disp_offset
            disp[disp > 0.96] -= 1.0
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
          output_flags):

    aa, aa_inv, nat, x_frac0 = read_OpenMX_input(dat_init)

    scale_disp, scale_force, scale_energy = get_unit_conversion_factor(
        str_unit)

    print_disp, print_force, print_energy, _ = output_flags

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
