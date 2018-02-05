#!/usr/bin/python
#
# extract_openMX_mod.py
#
# Simple script to extract atomic displacements, atomic forces, and
# energies from OpenMX output files.
# 
# Copyright (c) 2017 Yuto Tanaka
#
# Original script
# Copyright (c) 2014, 2015, 2016 Terumasa Tadano
# 
"""
This python script extracts atomic displacements, atomics forces,
and energies. 

---unit of data in nout files---
 * coordinatea  : angstrom
 * forces       : Hartree/Bohr
 * total energy : Hartree 
"""

import numpy as np


def get_original_OpenMX_mod(file_original):
    
    search_target1 = "Atoms.Number"
    search_target2 = "Atoms.SpeciesAndCoordinates.Unit"
    search_target3 = "<Atoms.SpeciesAndCoordinates"
    search_target4 = "<Atoms.UnitVectors"
    #open original file
    f = open(file_original, 'r')

    #set initial patameters
    nat = 0
    lavec_flag = 0
    lavec_row = 0
    lavec = np.zeros([3, 3])

    coord_flag = 0
    coord_row = 0

    #read oroginal file and pull out some infomations
    for line in f:
        ss = line.strip().split()
        #number of atoms
        if len(ss) > 0 and ss[0] == search_target1:
            nat = int(ss[1])
    
        #atomic coordinates
        if coord_flag == 1:
            for j in range(3):
                x_frac0[coord_row][j] = float(ss[j+2])

            coord_row += 1              
            if coord_row == nat:
                coord_flag = 0

        #latice vector
        if lavec_flag == 1:
            for i in range(3):
                lavec[lavec_row][i] = float(ss[i])
            lavec_row += 1
            if lavec_row == 3:
                lavec_flag = 0

        # unit of atomic coordinates
        if len(ss) > 0 and ss[0] == search_target2:
            coord_unit = ss[1]
        
        if len(ss) > 0 and ss[0] == search_target3:
            coord_flag = 1
            # initialize x_frac0 array
            x_frac0 = np.zeros([nat, 3])
 
        if len(ss) > 0 and ss[0] == search_target4:
            lavec_flag = 1

        if np.linalg.norm(lavec) > 0 and lavec_flag == 0:
            break


    #errors
    if nat == 0:
        print "Could not read dat file properly."
        exit(1)

      
    lavec_inv = (np.linalg.inv(lavec)).T
    #convert to frac
    if coord_unit == "ang" or coord_unit == "ANG" or coord_unit == "Ang":
        for i in range(nat):
            x_frac0[i] = np.dot(lavec_inv , x_frac0[i])
            
    f.close()

    return lavec, lavec_inv, nat, x_frac0



"""displacements"""

def get_coordinates_OpenMX(md_file, nat, lavec, conv):
    
    x = np.zeros([nat, 3], dtype = np.float64)

    f = open(md_file, 'r') 

    line_atom = f.readline()
    atom_c = int(line_atom.strip().split()[0])
    if atom_c != nat:
        print "The number of atom does not match."

    line = f.readline()

    atom_count = 0
    for line in f:
        ss = line.strip().split()
        for i in range(3):
            x[atom_count][i] = float(ss[i+1])
            
        # convert unit ang to frac
        x[atom_count] = np.dot(conv, x[atom_count]) 
        atom_count += 1

        if atom_count == nat:
            break

    f.close()
   
    return x


def print_displacements_OpenMX(md_files,
                        lavec, lavec_inv, nat, x0,
                        require_conversion,
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
        x0_offset = get_coordinates_OpenMX(file_offset, nat)
        try:
            x0_offset = np.reshape(x0_offset, (nat, 3))
        except:
            print"File %s contains too many position entries" % file_offset 
        disp_offset = x0_offset - x0

    for search_target in md_files:

        x = get_coordinates_OpenMX(search_target, nat, lavec, conv)
        #ndata = len(x) / (3 * nat)
        ndata = 1
        #x = np.reshape(x, (1, nat, 3))
       
        for idata in range(ndata):
            #disp = x[idata, :, :] - x0 - disp_offset
            disp = x - x0 - disp_offset
            disp[disp > 0.96] -= 1.0
            #disp = np.dot(vec_refold(disp), conv_inv)
            for i in range(nat):
                disp[i] = np.dot(conv_inv, disp[i])

            disp[np.absolute(disp) < 1e-5] = 0.0
            if require_conversion:
                disp *= conversion_factor

            for i in range(nat):
                print "%15.7F %15.7F %15.7F" % (disp[i][0],
                                                disp[i][1],
                                                disp[i][2]) 


"""atomic forces"""

def get_atomicforces_OpenMX(md_file, nat):

    force = np.zeros([nat, 3], dtype = np.float64)

    f = open(md_file, 'r')
    
    line_atom = f.readline()
    atom_c = int(line_atom.strip().split()[0])
    if atom_c != nat:
        print "The number of atom does not match." 

    line = f.readline()

    atom_count = 0
    for line in f:
        ss = line.strip().split()
        for i in range(3):
            force[atom_count][i] = float(ss[i+4])
        atom_count += 1

        if atom_count == nat:
            break

    f.close()

    return force


def print_atomicforces_OpenMX(md_files,
                            nat,
                            require_conversion,
                            conversion_factor,
                            file_offset):
  
    if file_offset is None:
        force_offset = np.zeros((nat, 3))
    else:
        data0 = get_atomicforces_OpenMX(file_offset, nat)
        try:
            force_offset = np.reshape(data0, (nat, 3))
        except:
            print "File %s contains too many force entries" % file_offset

    for search_target in md_files:
        data = get_atomicforces_OpenMX(search_target, nat)
        #ndata = len(data) / (3 * nat)
        ndata = 1
        #data = np.reshape(data, (ndata, nat, 3))

        for idata in range(ndata):
            #f = data[idata, :, :] - force_offset
            f = data - force_offset

            if require_conversion:
                f *= conversion_factor

            for i in range(nat):
                print "%15.8E %15.8E %15.8E" % (f[i][0],
                                                f[i][1],
                                                f[i][2]) 


"""total enegy"""
def get_energies_OpenMX(md_file, nat):

    target = "time="
    etot = []

    f = open(md_file, 'r')
    for line in f:
        ss = line.strip().split()
        if len(ss) > 0 and ss[0] == target:
            etot.extend([ss[4]])
            break
        else:
            continue

    if len(etot) == 0:
        print "Total energy not found." 
        exit(1)

    return np.array(etot, dtype=np.float)


# Other functions

def refold(x):
    if x >= 0.5:
        return x - 1.0
    elif x < -0.5:
        return x + 1.0
    else:
        return x


