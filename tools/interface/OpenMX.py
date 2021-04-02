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
                    nat = int(line.strip().split()[1])

                elif search_target[1] in line.lower():
                    coord_unit = line.strip().split()[1].lower()

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

            for line in lines[ipos_lavec:ipos_lavec + 3]:
                lavec.append([float(t) for t in line.strip().split()])

            common_settings.append(lines[:ipos_coord])
            common_settings.append(lines[fpos_coord:])

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

    def parse(self, initial_dat, out_files, out_file_offset, str_unit,
              output_flags, filter_emin=None, filter_emax=None):

        if not self._initial_structure_loaded:
            self.load_initial_structure(initial_dat)

        self._set_unit_conversion_factor(str_unit)
        self._set_output_flags(output_flags)

        if self._print_disp or self._print_force:
            self._print_displacements_and_forces(out_files,
                                                 out_file_offset,
                                                 filter_emin,
                                                 filter_emax)
        elif self._print_energy:
            self._print_energies(out_files, out_file_offset)

    def _generate_input(self, header, disp):

        filename = self._prefix + str(self._counter).zfill(self._nzerofills) + ".dat"

        with open(filename, 'w') as f:
            for line in self._common_settings[0]:
                if "atoms.speciesandcoordinates.unit" in line.lower():
                    f.write("Atoms.SpeciesAndCoordinates.Unit frac\n")
                else:
                    f.write("%s\n" % line)

            for i in range(self._nat):
                f.write("%4d %3s" % (i + 1, self._kd[i]))
                for j in range(3):
                    f.write("%20.16f" % (self._x_fractional[i, j] + disp[i, j]))
                for j in range(2):
                    f.write("%6.2f" % (self._initial_charges[i, j]))
                f.write('\n')

            for line in self._common_settings[1]:
                f.write("%s\n" % line)

        self._counter += 1

    def _print_displacements_and_forces(self, out_files,
                                        file_offset, filter_emin, filter_emax):

        # vec_refold = np.vectorize(refold)
        lavec_transpose = self._lattice_vector.transpose()

        x0 = np.round(self._x_fractional, 8)

        if file_offset is None:
            disp_offset = np.zeros((self._nat, 3))
            force_offset = np.zeros((self._nat, 3))
            epot_offset = 0.0

        else:
            x0_offset, force_offset = self._get_coordinate_and_force_outfile(file_offset)
            try:
                x0_offset = np.reshape(x0_offset, (self._nat, 3))
            except:
                raise RuntimeError("File %s contains too many position entries" % file_offset)

            disp_offset = x0_offset - x0
            try:
                force_offset = np.reshape(force_offset, (self._nat, 3))
            except:
                raise RuntimeError("File %s contains too many force entries" % file_offset)

            epot_offset = self._get_energies_outfile(file_offset)

        for search_target in out_files:

            x, force = self._get_coordinate_and_force_outfile(search_target)
            epot = self._get_energies_outfile(search_target)
            epot -= epot_offset
            epot *= self._RYDBERG_TO_EV * 2.0

            ndata = 1

            for idata in range(ndata):

                if filter_emin is not None:
                    if filter_emin > epot[idata]:
                        continue

                if filter_emax is not None:
                    if filter_emax < epot[idata]:
                        continue

                if self._print_disp:
                    disp = x - x0 - disp_offset
                    disp[disp > 0.96] -= 1.0
                    for i in range(self._nat):
                        disp[i] = np.dot(disp[i], lavec_transpose)

                    disp[np.absolute(disp) < 1e-5] = 0.0
                    disp *= self._disp_conversion_factor

                if self._print_force:
                    f = force - force_offset
                    f *= self._force_conversion_factor

                print("# Filename: %s, Snapshot: %d, E_pot (eV): %s" %
                      (search_target, idata + 1, epot[idata]))

                if self._print_disp and self._print_force:
                    for i in range(self._nat):
                        print("%15.7F %15.7F %15.7F %20.8E %15.8E %15.8E" % (disp[i, 0],
                                                                             disp[i, 1],
                                                                             disp[i, 2],
                                                                             f[i, 0],
                                                                             f[i, 1],
                                                                             f[i, 2]))
                elif self._print_disp:
                    for i in range(self._nat):
                        print("%15.7F %15.7F %15.7F" % (disp[i, 0],
                                                        disp[i, 1],
                                                        disp[i, 2]))
                elif self._print_force:
                    for i in range(self._nat):
                        print("%15.8E %15.8E %15.8E" % (f[i, 0],
                                                        f[i, 1],
                                                        f[i, 2]))

    def _print_energies(self, out_files, file_offset):

        if file_offset is None:
            etot_offset = 0.0
        else:
            etot_offset = self._get_energies_outfile(file_offset)

        print("# Etot")
        for search_target in out_files:

            etot = self._get_energies_outfile(search_target)

            for idata in range(len(etot)):
                val = etot[idata] - etot_offset
                val *= self._energy_conversion_factor
                print("%19.11E" % val)

    def _set_number_of_zerofill(self, npattern):

        nzero = 1

        while True:
            npattern //= 10
            if npattern == 0:
                break
            nzero += 1

        self._nzerofills = nzero

    def _set_unit_conversion_factor(self, str_unit):

        if str_unit == "ev":
            disp_conv_factor = 1.0
            energy_conv_factor = 2.0 * self._RYDBERG_TO_EV
            force_conv_factor = energy_conv_factor / self._BOHR_TO_ANGSTROM

        elif str_unit == "rydberg":
            disp_conv_factor = 1.0 / self._BOHR_TO_ANGSTROM
            energy_conv_factor = 2.0
            force_conv_factor = 2.0

        elif str_unit == "hartree":
            disp_conv_factor = 1.0 / self._BOHR_TO_ANGSTROM
            energy_conv_factor = 1.0
            force_conv_factor = 1.0

        else:
            raise RuntimeError("This cannot happen")

        self._disp_conversion_factor = disp_conv_factor
        self._force_conversion_factor = force_conv_factor
        self._energy_conversion_factor = energy_conv_factor

    def _set_output_flags(self, output_flags):
        self._print_disp, self._print_force, \
        self._print_energy, self._print_born = output_flags

    @property
    def nat(self):
        return self._nat

    @property
    def lattice_vector(self):
        return self._lattice_vector

    @property
    def inverse_lattice_vector(self):
        return self._inverse_lattice_vector

    @property
    def atomic_kinds(self):
        return self._kd

    @property
    def x_fractional(self):
        return self._x_fractional

    def _get_coordinate_and_force_outfile(self, out_file):
        """
        Return fractional coordinates and atomic forces in units of Hartree/Bohr
        """
        search_flag = "<coordinates.forces"
        f = open(out_file, 'r')
        line = f.readline()
        found_tag = False

        x = np.zeros((self._nat, 3))
        force = np.zeros((self._nat, 3))

        while line:
            if search_flag in line:
                found_tag = True
                f.readline()  # skip one line
                for i in range(self._nat):
                    line = f.readline()
                    x[i][:] = [float(t) for t in line.rstrip().split()[2:5]]
                    force[i][:] = [float(t) for t in line.rstrip().split()[5:]]
                break
            line = f.readline()

        if not found_tag:
            raise RuntimeError("%s tag not found in %s" % (search_flag, out_file))

        x = np.array(x)
        for i in range(self._nat):
            x[i, :] = np.dot(x[i, :], self._inverse_lattice_vector.transpose())

        return x, np.array(force)

    @staticmethod
    def _get_energies_outfile(out_file):

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
            raise RuntimeError("Total energy not found.")

        return np.array(etot, dtype=np.float)
