#!/usr/bin/env python

import sys
import numpy as np
import optparse
import math

parser = optparse.OptionParser()
parser.add_option('-t', '--temp')
parser.add_option('-m', '--mode')
options, args = parser.parse_args()

file_result = args[0]

fs = open(file_result, 'r')


def locate_tag(tag, fs):

    fs.seek(0)
    key_exist = False

    while True:
        line = fs.readline()

        if not line:
            break

        if line.rstrip() == tag:
            key_exist = True
            break

    if not key_exist:
        sys.exit("key " + tag + " not found")
        

if __name__ == '__main__':

    locate_tag("#SYSTEM", fs)
    data = fs.readline().split()
    nat, nkd = [int(data[i]) for i in range(len(data))]
    ns = 3 * nat

    locate_tag("#TEMPERATURE", fs)
    data = fs.readline().split()
    tmin, tmax, dt = [int(data[i]) for i in range(len(data))]
    nt = (tmax - tmin) / dt
    temp = np.zeros(nt)
    for i in range(nt):
        temp[i] = tmin + dt * float(i)

    locate_tag("#KPOINT", fs)
    data = fs.readline().split()
    nkx, nky, nkz = [int(data[i]) for i in range(len(data))]
    nk = int(fs.readline())

    nks = nk * ns

    locate_tag("##Phonon Frequency", fs)
    fs.readline()

    omega = np.zeros((nk, ns))
    tau = np.zeros((nk, ns, nt))
    vel_rep = np.zeros((nk, ns, 3))

    for i in range(nk):
        for j in range(ns):
            data = fs.readline().split()
            omega[i,j] = float(data[2])

    if options.temp == None:
        if options.mode == None:
            sys.exit("Not supported yet")
        else:
            print "# mode number", options.mode, "will be presented"
            target_k = int(options.mode) / ns
            target_s = int(options.mode) % ns
            calc = "damping"
    else:
        print "# Phonon properties at Temperature ", options.temp, "K"
        
        if (float(options.temp) - tmin) % dt != 0:
            sys.exit("No information is found at the given temperature")
        else:
            itemp = int((float(options.temp) - tmin) / dt)
        
        if options.mode == None:
            print "# All mode will be printed"
            calc = "tau_all"
        else:
            target_k = int(options.mode) / ns
            target_s = int(options.mode) % ns
            calc = "tau_mode"
            
    locate_tag("##Phonon Relaxation Time", fs)

    # Read Relaxation Times
    
    for i in range(nk):
        for j in range(ns):

            fs.readline()
            fs.readline()
            
            nk_equiv = int(fs.readline())
            
            # Read velocity

            for k in range(nk_equiv):
                if k == 0:
                    data = fs.readline().split()
                    vel_rep[i,j,0], vel_rep[i,j,1], vel_rep[i,j,2] = [float(data[l]) for l in range(3)]
                else:
                    fs.readline()

            for k in range(nt):
                tau[i, j, k] = float(fs.readline())
        
            fs.readline()

    
    # Write phonon info

    if calc == "damping":
        print "# Mode numbers", target_k + 1, target_s + 1
        print "# Temperature [k], Relaxation Time [ps], MFP [nm]"

        vel_norm = math.sqrt(vel_rep[target_k, target_s, 0]**2 + vel_rep[target_k, target_s, 1]**2 + vel_rep[target_k, target_s, 2]**2)
        
        for i in range(nt):
            print "%9.3f %15.7f %15.7f" % (temp[i], tau[target_k, target_s, i], tau[target_k, target_s, i]*vel_norm*0.001)

    elif calc == "tau_mode":
        print "# Mode numbers", target_k + 1, target_s + 1
        print "# Frequency = ", omega[target_k, target_s], " cm^{-1}"
        print "# Temperature [k], Relaxation Time [ps], velocity [m/s],  MFP [nm]"


        vel_norm = math.sqrt(vel_rep[target_k, target_s, 0]**2 + vel_rep[target_k, target_s, 1]**2 + vel_rep[target_k, target_s, 2]**2)
        print "%9.3f %15.7f %15.7f %15.7f" % (temp[itemp], tau[target_k, target_s, itemp], vel_norm, tau[target_k, target_s, itemp]*vel_norm*0.001)
        
    elif calc == "tau_all":
        print "# ik, is, Frequency [cm^{-1}], Relaxation Time [ps], MFP [nm]"

        for i in range(nk):
            for j in range(ns):
                vel_norm = math.sqrt(vel_rep[i,j,0]**2 + vel_rep[i,j,1]**2 + vel_rep[i,j,2]**2)
                print "%5i %5i %9.3f %15.7f %15.7f %15.7f" % (i + 1, j + 1, omega[i,j], tau[i,j,itemp], vel_norm, tau[i,j,itemp]*vel_norm*0.001)
                
        

