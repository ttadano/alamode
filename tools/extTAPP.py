#!/usr/bin/env python
#
# extTAPP.py
#
# Simple script to extract atomic displacements, atomic forces, and
# energies from xTAPP *.str files
#
# Copyright (c) 2014 Terumasa Tadano
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory 
# or http://opensource.org/licenses/mit-license.php for information.
#
"""
This python script extracts atomic displacements, atomics forces,
and energies from xTAPP *.str files.
"""

import numpy as np
import sys, os, math
import optparse

usage = "usage: %prog [options] *.str"
parser = optparse.OptionParser(usage=usage)
parser.add_option('--unit', action="store", type="string", dest="unitname", default="Rydberg",
    help="print atomic displacements and forces in units of UNIT. Available options are eV, Rydberg (default), and Hartree.")
parser.add_option('--get', help="specify which quantity to extract. Available options are disp, force and energy.")
parser.add_option('--cg', help="Original CG file with equilibrium atomic positions (default: None)")


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
        print "main data entry not found"
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
            a = float(entrylist[i+2])

        if "lattice_list" in entrylist[i]:
            for j in range(9):
                lavec_list.append(entrylist[i+j+2])

        if "number_element" in entrylist[i]:
            nkd = int(entrylist[i+2])

        if "number_atom" in entrylist[i]:
            nat = int(entrylist[i+2])
    
    if a == 0.0:
        print "Couldn't read lattice_factor"
        exit(1)
    if nkd == 0:
        print "Couldn't read number_element"
        exit(1)
    if nat == 0:
        print "Couldn't read number_atom"
        exit(1)
    if len(lavec_list) != 9:
        print "Couldn't read lattice_list"
        exit(1)

    lavec = np.zeros((3,3))

    for i in range(3):
        for j in range(3):
            lavec[j][i] = a * float(lavec_list[3 * i + j])

    return lavec, nat, nkd



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
        print "atom data entry not found"
        exit(1)

    x_out = np.zeros((nat_in, 3), dtype=float)
    kd_out = np.zeros(nat_in, dtype=int)

    for i in range(nat_in):
        list_tmp = list_atom[i + nkd_in + 1].rstrip().split()
        kd_out[i] = int(list_tmp[0])
        for j in range(3):
            x_out[i][j] = float(list_tmp[j+1])

    return x_out, kd_out

def read_CG(file_in):
    
    lavec, nat, nkd = read_tappinput(file_in)
    x0, kd = read_atomdata(file_in, nat, nkd)

    lavec = np.matrix(lavec)
    lavec_inv = np.array(lavec.I)
    lavec = np.array(lavec)

    return lavec, nat, x0


def refold(x):
    if x >= 0.5:
        return x - 1.0
    elif x < -0.5:
        return x + 1.0
    else:
        return x


def print_displacements(str_files, lavec, nat, x0, require_conversion, conversion_factor):


    x0 = np.round(x0, 8)
    x = np.zeros((nat, 3))

    lavec_transpose = lavec.T

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
            print "atom_position tag not found in %s" % search_target
            exit(1)


        disp = x - x0
        for i in range(nat):
            disp[i,:] = [refold(disp[i,j]) for j in range(3)]

        disp = np.dot(disp, lavec_transpose)

        if require_conversion:
            disp *= conversion_factor

        for i in range(nat):
            print "%15.7F %15.7F %15.7F" % (disp[i][0], disp[i][1], disp[i][2])


def print_atomicforces(str_files, nat, require_conversion, conversion_factor):

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
            print "force tag not found in %s" % search_target
            exit(1)


        if require_conversion:
            force *= conversion_factor

        for i in range(nat):
            print "%19.11E %19.11E %19.11E" % (force[i][0], force[i][1], force[i][2])
    

def print_energies(str_files, require_conversion, conversion_factor):

    for search_target in str_files:

        found_tag = False

        with open(search_target) as openfileobject:
            for line in openfileobject:
                if "total_energy" in line:
                    energy_str = line.rstrip().split()[2]
                    etot =  float(energy_str[:-1])

                    if require_conversion:
                        etot *= conversion_factor

                    print "%19.11E" % etot

                    found_tag = True
                    break
        
        if not found_tag:
            print "total_energy tag not found in %s" % search_target
            exit(1)


if __name__ == "__main__":

    options, args = parser.parse_args()
    file_str = args[0:]

    if len(file_str) == 0:
        print "Usage: extTAPP.py [options] *.str"
        print 
        print "For details of available options, please type\n$ python extTAPP.py -h"
        exit(1)

    if options.cg == None:
        print "Original CG file is not specified by -i option"
        exit(1)
    else:
        file_cg = options.cg

    Bohr_radius = 0.52917721092 #Angstrom
    Rydberg_to_eV = 13.60569253

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
        print "Invalid options for --unit"
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

    aa, nat, x_frac0 = read_CG(file_cg)

    if print_disp:
        print_displacements(file_str, aa, nat, x_frac0, convert_unit, disp_conv_factor)
    elif print_force:
        print_atomicforces(file_str, nat, convert_unit, force_conv_factor)
    elif print_energy:
        print_energies(file_str, convert_unit, energy_conv_factor)




