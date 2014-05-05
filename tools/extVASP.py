#!/usr/bin/env python
#
# extVASP.py
#
# Simple script to extract atomic displacements, atomic forces, and
# energies from vasprun.xml files
#
# Copyright (c) 2014 Terumasa Tadano
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory 
# or http://opensource.org/licenses/mit-license.php for information.
#
"""
This python script extracts atomic displacements, atomics forces,
and energies from vasprun.xml files.
"""

import numpy as np
import sys, os, math
import optparse

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

usage = "usage: %prog [options]  */vasprun.xml"
parser = optparse.OptionParser(usage=usage)
parser.add_option('-u', '--unit', action="store", type="string", dest="unitname", default="Rydberg",
    help="print atomic displacements and forces in units of UNIT. Available options are original, Rydberg (default), and Hartree.")
parser.add_option('-e', '--extract', help="specify which quantity to extract. Available options are disp, force and energy.")
parser.add_option('-i', '--input', help="Original POSCAR file with equilibrium atomic positions (default: POSCAR)")


def read_POSCAR(file_in):
    file_pos = open(file_in, 'r')
    
    str_tmp = file_pos.readline()
    a = float(file_pos.readline().rstrip())
    lavec = np.zeros((3,3))

    for i in range(3):
        arr = file_pos.readline().rstrip().split()
        if len(arr) != 3:
            print "Could not read POSCAR properly"
            exit(1)

        for j in range(3):
            lavec[i,j] = a * float(arr[j])

    lavec = np.matrix(lavec).transpose()
    invlavec = lavec.I

    elements = file_pos.readline().rstrip().split()
    nat_elem = [int(tmp) for tmp in file_pos.readline().rstrip().split()]

    nat = np.sum(nat_elem)
    basis = file_pos.readline().rstrip()
    x = np.zeros((nat, 3))

    for i in range(nat):
        arr = file_pos.readline().rstrip().split()
        for j in range(3):
            x[i,j] = float(arr[j])

    if basis == "Direct" or basis == "direct" or basis == "D" or basis == "d":
        xf = np.matrix(x)
    else:
        xf = np.matrix(x)
        for i in range(nat):
            xf[i,:] = xf[i,:] * invlavec.transpose()

    file_pos.close()
    return lavec, invlavec, elements, nat_elem, xf


def refold(x):
    if x >= 0.5:
        return x - 1.0
    elif x < -0.5:
        return x + 1.0
    else:
        return x


def print_displacements(xml_files, lavec, nat, x0, require_conversion, conversion_factor):

    for i in range(nat):
        for j in range(3):
            x0[i, j] = round(x0[i, j], 8)


    x = np.zeros((nat, 3))

    lavec_transpose = lavec.transpose()

    for search_target in xml_files:

        xml = etree.parse(search_target)
        root = xml.getroot()

        for elems in root.findall('calculation/structure/varray'):
            str_coord = []
            for elems2 in elems.findall('v'):
                str_coord.append(elems2.text)
                    
            n = len(str_coord)

            for i in range(n):
                x[i, :] = [float(t) for t in str_coord[i].split()]
            
            x = np.matrix(x)
            disp = x - x0

            for i in range(n):
                disp[i,:]= [refold(disp[i,j]) for j in range(3)]

            disp = disp * lavec_transpose

            if require_conversion:
                for i in range(n):
                    disp[i,:] *= conversion_factor

            for i in range(n):
                print "%15.7F %15.7F %15.7F" % (disp[i,0], disp[i,1], disp[i,2])


def print_atomicforces(xml_files, nat, require_conversion, conversion_factor):

    f = np.zeros((nat, 3))
    for search_target in xml_files:

        xml = etree.parse(search_target)
        root = xml.getroot()

        for elems in root.findall('calculation/varray'):
            str_force = []
            if elems.get('name') == "forces":
                for elems2 in elems.findall('v'):
                    str_force.append(elems2.text)
                    
            n = len(str_force)

            for i in range(n):
                f[i, :] = [float(t) for t in str_force[i].split()]

            if require_conversion:
                for i in range(n):
                    f[i, :] *= conversion_factor

            for i in range(n):
                 print "%15.8E %15.8E %15.8E" % (f[i,0], f[i,1], f[i,2])

def print_energies(xml_files, require_conversion, conversion_factor):

    print "# Etot, Ekin"

    for search_target in xml_files:

        xml = etree.parse(search_target)
        root = xml.getroot()

        for elems in root.findall('calculation/energy'):
            etot = 'N/A'
            ekin = 'N/A'
            
            for elems2 in elems.findall('i'):
                if elems2.get('name') == "e_fr_energy":
                    etot = elems2.text
                if elems2.get('name') == "kinetic":
                    ekin = elems2.text
            
            if not require_conversion:
                print "%s %s" % (etot, ekin)
            else:
                if etot != 'N/A':
                    val_etot = float(etot) * conversion_factor
                    print "%15.8E" % val_etot,
                else:
                    print "%s" % etot,
                
                if ekin != 'N/A':
                    val_ekin = float(ekin) * conversion_factor
                    print "%15.8E" % val_ekin
                else:
                    print "%s" % ekin


if __name__ == "__main__":

    options, args = parser.parse_args()
    file_xml = args[0:]

    if len(file_xml) == 0:
        print "Usage: extVASP.py [options] */vasprun.xml"
        print 
        print "For details of available options, please type\n$ python extVASP.py -h"
        exit(1)

    if options.input == None:
        file_poscar = "POSCAR"
    else:
        file_poscar = options.input

    Bohr_radius = 0.52917721092 #Angstrom
    Rydberg_to_eV = 13.60569253

    if options.unitname == "original":
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
        print "Invalid options for --unit"
        exit(1)

    force_conv_factor = energy_conv_factor / disp_conv_factor

    print_disp = False
    print_force = False
    print_energy = False

    if options.extract == "disp":
        print_disp = True
    elif options.extract == "force":
        print_force = True
    elif options.extract == "energy":
        print_energy = True
    else:
        print "Please specify which quantity to extract by the --extract option."
        exit(1)

    aa, aa_inv, elems, nats, x_frac0 = read_POSCAR(file_poscar)
    #print "Total number of atoms: ", np.sum(nats)



    if print_disp:
        print_displacements(file_xml, aa, np.sum(nats), x_frac0, convert_unit, disp_conv_factor)
    elif print_force:
        print_atomicforces(file_xml, np.sum(nats), convert_unit, force_conv_factor)
    elif print_energy:
        print_energies(file_xml, convert_unit, energy_conv_factor)


  

  
