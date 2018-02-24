#
# VASP.py
#
# Interface to VASP (https://www.vasp.at)
#
# Copyright (c) 2014 Terumasa Tadano
#
# This file is distributed under the terms of the MIT license.
# Please see the file 'LICENCE.txt' in the root directory
# or http://opensource.org/licenses/mit-license.php for information.
#
from __future__ import print_function
import numpy as np

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

def read_POSCAR(file_in):

    file_pos = open(file_in, 'r')

    file_pos.readline()
    a = float(file_pos.readline().rstrip())
    lavec = np.zeros((3, 3))

    for i in range(3):
        arr = file_pos.readline().rstrip().split()
        if len(arr) != 3:
            print("Could not read POSCAR properly")
            exit(1)

        for j in range(3):
            lavec[i, j] = a * float(arr[j])

    lavec = lavec.transpose()
    invlavec = np.linalg.inv(lavec)

    elements = file_pos.readline().rstrip().split()

    if elements[0].isdigit():
        nat_elem = [int(tmp) for tmp in elements]
        elements = []

    else:
        nat_elem = [int(tmp) for tmp in file_pos.readline().rstrip().split()]

    nat = np.sum(nat_elem)
    basis = file_pos.readline().rstrip()
    x = np.zeros((nat, 3))

    for i in range(nat):
        arr = file_pos.readline().rstrip().split()
        for j in range(3):
            x[i][j] = float(arr[j])

    if basis == "Direct" or basis == "direct" or basis == "D" or basis == "d":
        xf = x
    else:
        xf = np.dot(x, invlavec)

    file_pos.close()

    return lavec, invlavec, elements, nat_elem, xf


def write_POSCAR(prefix, counter, header, nzerofills,
                 lavec, elems, nat, disp, coord):

    filename = prefix + str(counter).zfill(nzerofills) + ".POSCAR"
    f = open(filename, 'w')
    f.write("%s\n" % header)
    f.write("%s\n" % "1.0")

    for i in range(3):
        f.write("%20.15f %20.15f %20.15f\n" % (lavec[0][i],
                                               lavec[1][i],
                                               lavec[2][i]))

    for i in range(len(elems)):
        f.write("%s " % elems[i])
    if len(elems) > 0:
        f.write("\n")

    for i in range(len(nat)):
        f.write("%d " % nat[i])
    f.write("\n")

    f.write("Direct\n")

    for i in range(len(disp)):
        for j in range(3):
            f.write("%20.15f" % (coord[i][j] + disp[i][j]))
        f.write("\n")
    f.close()




def get_coordinate_VASP(xml_file, nat):

    x = []

    try:
        xml = etree.parse(xml_file)
        root = xml.getroot()

        for elems in root.findall('calculation/structure/varray'):
            str_coord = [elems2.text for elems2 in elems.findall('v')]
            n = len(str_coord)

            for i in range(n):
                x.extend([t for t in str_coord[i].split()])

        return np.array(x, dtype=np.float)

    except:
        print("Error in reading atomic positions from the XML file: %s" % xml_file)


def print_displacements_VASP(xml_files,
                             lavec, nat, x0,
                             require_conversion,
                             conversion_factor,
                             file_offset):

    x0 = np.round(x0, 8)
    lavec_transpose = lavec.transpose()
    vec_refold = np.vectorize(refold)

    if file_offset is None:
        disp_offset = np.zeros((nat, 3))
    else:
        x0_offset = get_coordinate_VASP(file_offset, nat)
        try:
            x0_offset = np.reshape(x0_offset, (nat, 3))
        except:
            print("File %s contains too many position entries" % file_offset)
        disp_offset = x0_offset - x0

    for search_target in xml_files:

        x = get_coordinate_VASP(search_target, nat)
        ndata = len(x) // (3 * nat)
        x = np.reshape(x, (ndata, nat, 3))

        for idata in range(ndata):
            disp = x[idata, :, :] - x0 - disp_offset
            disp = np.dot(vec_refold(disp), lavec_transpose)

            if require_conversion:
                disp *= conversion_factor

            for i in range(nat):
                print("%15.7F %15.7F %15.7F" % (disp[i, 0],
                                                disp[i, 1],
                                                disp[i, 2]))


def get_atomicforces_VASP(xml_file):

    f = []

    try:
        xml = etree.parse(xml_file)
        root = xml.getroot()

        for elems in root.findall('calculation/varray'):
            if elems.get('name') == "forces":
                str_force = [elems2.text for elems2 in elems.findall('v')]

                for i in range(len(str_force)):
                    f.extend([t for t in str_force[i].split()])
                
        return np.array(f, dtype=np.float)

    except:
        print("Error in reading atomic forces from the XML file: %s" % xml_file)


def print_atomicforces_VASP(xml_files,
                            nat,
                            require_conversion,
                            conversion_factor,
                            file_offset):

    if file_offset is None:
        force_offset = np.zeros((nat, 3))
    else:
        data0 = get_atomicforces_VASP(file_offset)
        try:
            force_offset = np.reshape(data0, (nat, 3))
        except:
            print("File %s contains too many force entries" % file_offset)

    for search_target in xml_files:

        data = get_atomicforces_VASP(search_target)
        ndata = len(data) // (3 * nat)
        data = np.reshape(data, (ndata, nat, 3))

        for idata in range(ndata):
            f = data[idata, :, :] - force_offset

            if require_conversion:
                f *= conversion_factor

            for i in range(nat):
                print("%15.8E %15.8E %15.8E" % (f[i][0],
                                                f[i][1],
                                                f[i][2]))


def get_energies_VASP(xml_file):

    etot_array = []
    ekin_array = []

    try:
        xml = etree.parse(xml_file)
        root = xml.getroot()

        for elems in root.findall('calculation/energy'):
            etot = 'N/A'
            ekin = 'N/A'

            for elems2 in elems.findall('i'):
                if elems2.get('name') == "e_fr_energy":
                    etot = elems2.text
                if elems2.get('name') == "kinetic":
                    ekin = elems2.text

            etot_array.append(etot)
            ekin_array.append(ekin)

        return etot_array, ekin_array
    except:
        print("Error in reading energies from the XML file: %s" % xml_file)


def print_energies_VASP(xml_files,
                        require_conversion,
                        conversion_factor,
                        file_offset):

    print("# Etot, Ekin")

    etot_offset = 0.0
    ekin_offset = 0.0

    if file_offset:
        etot, ekin = get_energies_VASP(file_offset)
        if len(etot) > 1 or len(ekin) > 1:
            print("File %s contains too many energy entries" % file_offset)
            exit(1)
        if etot[0] != 'N/A':
            etot_offset = float(etot[0])
        if ekin[0] != 'N/A':
            ekin_offset = float(ekin[0])

    for search_target in xml_files:

        etot, ekin = get_energies_VASP(search_target)

        for i in range(len(etot)):
            if etot[i] != 'N/A':
                val_etot = float(etot[i]) - etot_offset
                if require_conversion:
                    print("%15.8E" % (val_etot * conversion_factor), end=' ')
                else:
                    print("%15.8E" % val_etot, end=' ')
            else:
                print("%s" % etot[i], end=' ')

            if ekin[i] != 'N/A':
                val_ekin = float(ekin[i]) - ekin_offset
                if require_conversion:
                    print("%15.8E" % (val_ekin * conversion_factor))
                else:
                    print("%15.8E" % val_ekin)
            else:
                print("%s" % ekin[i])



def refold(x):
    if x >= 0.5:
        return x - 1.0
    elif x < -0.5:
        return x + 1.0
    else:
        return x

