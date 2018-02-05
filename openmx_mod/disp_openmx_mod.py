#!/usr/bin/env python
#
# displace.py
#
# Simple script to generate OpenMX input files 
# of given displacement patterns.
#
# Copyright (c) 2017 Yuto Tanaka
#
# oroginal script
# Copyright (c) 2014, 2015, 2016 Terumasa Tadano
#

"""
Input file generator for displaced configurations.
"""

import numpy as np
#from collections import OrderedDict

def read_OpenMX_input(file_in):

    search_target1 = "Atoms.Number"
    search_target2 = "<Atoms.UnitVectors"
    #open original file
    f = open(file_in, 'r')

    #set initial patameters
    nat = 0
    lavec_flag = 0
    lavec_row = 0
    lavec = np.zeros((3, 3))


    #read oroginal file and pull out some infomations
    for line in f:
        ss = line.strip().split()
        #number of atoms
        if len(ss) > 0 and ss[0] == search_target1:
            nat = int(ss[1])
 
        #latice vector
        if lavec_flag == 1:
            for i in range(3):
                lavec[lavec_row][i] = float(ss[i])
            lavec_row += 1
            if lavec_row == 3:
                lavec_flag = 0
 
        if len(ss) > 0 and ss[0] == search_target2:
            lavec_flag = 1

        if np.linalg.norm(lavec) > 0 and lavec_flag == 0:
            break
            
    f.close()

    #errors
    if nat == 0:
        print "Could not read dat file properly."
        exit(1)

    #calculate reciprocal vector
    lavec_T = lavec.transpose()
    invlavec = np.linalg.inv(lavec_T)

    return nat, lavec, invlavec


def write_OpenMX_input(prefix, counter, nzerofills, disp, lavec, file_in):
  
    search_target1 = "Atoms.Number"
    search_target2 = "<Atoms.SpeciesAndCoordinates"
    search_target3 = "Atoms.SpeciesAndCoordinates.Unit"
    search_target4 = "System.Name"

    str_ang = ["ANG", "ang", "Ang"]

    filename = prefix + str(counter).zfill(nzerofills) + ".dat"
    fout = open(filename, 'w')
    fin = open(file_in, 'r')

    nat = 0
    coord_flag = 0
    coord_row = 0

    conv = (np.linalg.inv(lavec)).T
    conv_inv = np.linalg.inv(conv)
    
    for i in range(nat):
        print np.dot(conv_inv, disp[i]) 
    
    disp[disp < 0] += 1

    for line in fin:
        ss = line.strip().split()
        #number of atoms
        if len(ss) > 0 and ss[0] == search_target1:
            nat = int(ss[1])
            x_frac = np.zeros((nat, 3))
            #coord = OrderedDict()
            coord = {}
            for i in range(nat):
                coord[i+1] = []
 
        #coordinates_unit
        if len(ss) > 0 and ss[0] == search_target3:
            coord_unit = ss[1]
        
        #coordinates
        if coord_flag == 1:
            coord_column = len(ss)
            for i in range(1, coord_column):
                if i > 1:
                    coord[int(ss[0])].append(float(ss[i]))
                else:
                    coord[int(ss[0])].append(ss[i])
       
            #convert to frac
            if coord_unit in str_ang:
                coord[coord_row+1] = np.dot(conv, coord[coord_row+1])

            # add displacement
            for j in range(1, 4):
                coord[coord_row+1][j] += disp[coord_row][j-1]
                coord[coord_row+1][j] = format(coord[coord_row+1][j],'20.16f')

            fout.write(str(coord_row+1) + " ")
            fout.write("  ".join(map(str, coord[coord_row+1])))
            fout.write("\n")
            coord_row += 1
            if coord_row == nat:
                coord_flag = 0

        elif len(ss) > 0 and ss[0] == search_target4:
            ss[1] = prefix + str(counter).zfill(nzerofills)
            fout.write("                      ".join(map(str, ss)))
            fout.write("\n")

        else:
            fout.write(line)


        if len(ss) > 0 and ss[0] == search_target2:
            coord_flag = 1


    fin.close()
    fout.close()


