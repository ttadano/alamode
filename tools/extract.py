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
import interface.VASP as vasp
import interface.QE as qe
import interface.xTAPP as xtapp
import interface.OpenMX as openmx
import interface.LAMMPS as lammps

usage = "usage: %prog [options] vasprun*.xml (or *.pw.out or *.str or *.md or *.out)"
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

parser.add_option('--OpenMX',
                  metavar='orig.dat',
                  help="OpenMX dat file with equilibrium atomic \
                        positions (default: None)")

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
                  options.LAMMPS is None,
                  options.OpenMX is None]

    if conditions.count(True) == len(conditions):
        print("Error : Either --VASP, --QE, --xTAPP, --LAMMPS, --OpenMX option must be given.")
        exit(1)

    elif len(conditions) - conditions.count(True) > 1:
        print("Error : --VASP, --QE, --xTAPP, --LAMMPS, and --OpenMX cannot be given simultaneously.")
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

        elif options.unitname == "Rydberg":
            convert_unit = True
            disp_conv_factor = 1.0 / Bohr_radius
            energy_conv_factor = 1.0 / Rydberg_to_eV

        elif options.unitname == "Hartree":
            convert_unit = True
            disp_conv_factor = 1.0 / Bohr_radius
            energy_conv_factor = 0.5 / Rydberg_to_eV
        
        else:
            print("Error : Invalid option for --unit")
            exit(1)

    
    elif options.OpenMX:
        code = "OpenMX"
        file_original = options.OpenMX

        if options.unitname == "eV":
            convert_unit = True
            disp_conv_factor = 1.0
            energy_conv_factor = 2.0 * Rydberg_to_eV
            force_conv_factor = energy_conv_factor 

        elif options.unitname == "Rydberg":
            convert_unit = True
            disp_conv_factor = 1.0 / Bohr_radius
            energy_conv_factor = 2.0
            force_conv_factor = 2.0

        elif options.unitname == "Hartree":
            convert_unit = True
            disp_conv_factor = 1.0 / Bohr_radius
            energy_conv_factor = 1.0
            force_conv_factor = 1.0

        else:
            print("Error : Invalid option for --unit")
            exit(1)

    if options.OpenMX is None:
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
        aa, aa_inv, elems, nats, x_frac0 = vasp.read_POSCAR(file_original)

    elif code == "QE":
        alat, aa, nat, x_frac0 = qe.read_original_QE_mod(file_original)

    elif code == "xTAPP":
        aa, nat, x_frac0 = xtapp.read_CG_mod(file_original)
    
    elif code == "LAMMPS":
        common_settings, nat, x_cart0, kd = lammps.read_lammps_structure(file_original)
    
    elif code == "OpenMX":
        aa, aa_inv, nat, x_frac0 = openmx.read_OpenMX_input(file_original)


    # Print data

    if print_disp:
        if code == "VASP":
            vasp.print_displacements_VASP(file_results, aa, np.sum(nats), x_frac0,
                                     convert_unit, disp_conv_factor, options.offset)

        elif code == "QE":
            qe.print_displacements_QE(file_results, alat, aa, nat, x_frac0,
                                   convert_unit, disp_conv_factor, options.offset)

        elif code == "xTAPP":
            xtapp.print_displacements_xTAPP(file_results, aa, nat, x_frac0,
                                      convert_unit, disp_conv_factor, options.offset)

        elif code == "LAMMPS":
            lammps.print_displacements_LAMMPS(file_results, nat, x_cart0,
                                       convert_unit, disp_conv_factor, options.offset)
    
        elif code == "OpenMX":
            openmx.print_displacements_OpenMX(file_results, aa, aa_inv, nat, x_frac0,
                                      convert_unit, disp_conv_factor, options.offset)

    elif print_force:
        if code == "VASP":
            vasp.print_atomicforces_VASP(file_results, np.sum(nats),
                                    convert_unit, force_conv_factor, options.offset)

        elif code == "QE":
            qe.print_atomicforces_QE(file_results, nat,
                                  convert_unit, force_conv_factor, options.offset)

        elif code == "xTAPP":
            xtapp.print_atomicforces_xTAPP(file_results, nat,
                                     convert_unit, force_conv_factor, options.offset)

        elif code == "LAMMPS":
            lammps.print_atomicforces_LAMMPS(file_results, nat, 
                                      convert_unit, force_conv_factor, options.offset)

        elif code == "OpenMX":
            openmx.print_atomicforces_OpenMX(file_results, nat,
                                     convert_unit, force_conv_factor, options.offset)

    elif print_energy:
        if code == "VASP":
            vasp.print_energies_VASP(file_results, convert_unit,
                                energy_conv_factor, options.offset)

        elif code == "QE":
            qe.print_energies_QE(file_results, convert_unit,
                              energy_conv_factor, options.offset)

        elif code == "xTAPP":
            xtapp.print_energies_xTAPP(file_results, convert_unit,
                                 energy_conv_factor, options.offset)

        elif code == "OpenMX":
            openmx.print_energies_OpenMX(file_results, convert_unit,
                                 energy_conv_factor, options.offset)
