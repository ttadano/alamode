#!/usr/bin/env python

import sys
import datetime
import numpy as np
import optparse
import math

# CONSTANTS
SPEED_C = 299792458.0
RYD     = 4.35974394e-18 / 2.0
KB      = 1.3806488e-23
BOHR_IN_AA = 0.52917721092

HZ_TO_KAYSER = 1.0e-2 / (2.0*math.pi*SPEED_C)
TIME_RY = 6.62606896e-34/ (2.0*math.pi*RYD)
T_TO_RYD = KB / RYD
RYD_TO_KAYSER = HZ_TO_KAYSER / TIME_RY
KAYSER_TO_RYD = 1.0 / RYD_TO_KAYSER
#

parser = optparse.OptionParser()
parser.add_option('-t', '--temp')
parser.add_option('-m', '--mode')
parser.add_option('-k', '--kpoint')
parser.add_option('-c', '--calc')
parser.add_option('-l', '--length')

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


def Cv(omega, T):
    if abs(omega) < 1.0e-12 or T == 0.0:
        return 0.0
    else:
        x = omega*KAYSER_TO_RYD / (T*T_TO_RYD)
        return KB * (x/(2.0 * math.sinh(0.5*x)))**2
        

if __name__ == '__main__':

    volume = float(args[1])

    print "#", datetime.datetime.today()
    print "# Input file : ", file_result
    print "#"

    # Parse phonon properties

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
    vel_sum = np.zeros((nk, ns, 3, 3))

    factor = 1.0e+18 / (BOHR_IN_AA**3 * float(nkx*nky*nkz) * volume)

    for i in range(nk):
        for j in range(ns):
            data = fs.readline().split()
            omega[i,j] = float(data[2])

    # Choose what to calculate

    if options.calc == None or options.calc == "kappa":
        calc = "kappa"
    elif options.calc == "kappa_size":
        calc = "kappa_size"
    elif options.calc == "tau":
        calc = "tau"
    else:
        sys.exit("Invalid --value option")
    
    if calc == "tau":
        if options.temp == None: # Damping function 
            if options.kpoint == None or options.mode == None:
                sys.exit("Invalid combination of --temp, --kpoint, and --mode options")
            else:
                if len(options.kpoint.split(':')) != 1:
                    sys.exit("Invalid usage of --kpoint for --calc=tau")
                if len(options.mode.split(':')) != 1:
                    sys.exit("Invalid usage of --mode for --calc=tau")

                target_k = int(options.kpoint) - 1
                target_s = int(options.mode) - 1
                calc = "damping"

        else: # Relaxation time (omega,tau)
            if (float(options.temp) - tmin) % dt != 0:
                sys.exit("No information is found at the given temperature")
            else:
                itemp = int((float(options.temp) - tmin) / dt)

            if options.kpoint == None:
                beg_k = 0
                end_k = nk
            else:
                if len(options.kpoint.split(':')) == 1:
                    beg_k = int(options.kpoint) - 1
                    end_k = beg_k + 1
                elif len(options.kpoint.split(':')) == 2:
                    arr = options.kpoint.split(':')
                    beg_k, end_k = int(arr[0]) - 1, int(arr[1])
                else:
                    sys.exit("Invalid usage of --kpoint for --calc=tau")

            if options.mode == None:
                beg_s = 0
                end_s = ns
            else:
                if len(options.mode.split(':')) == 1:
                    beg_s = int(options.mode) - 1
                    end_s = beg_s + 1
                elif len(options.mode.split(':')) == 2:
                    arr = options.mode.split(':')
                    beg_s, end_s = int(arr[0]) - 1, int(arr[1]) 
                else:
                    sys.exit("Invalid usage of --mode for --calc=tau")

    elif calc == "kappa_size":
        if options.temp == None:
            sys.exit("--temp should be given when --value = kappa_size")
        else:
            if (float(options.temp) - tmin) % dt != 0:
                sys.exit("No information is found at the given temperature")
            else:
                itemp = int((float(options.temp) - tmin) / dt)
        
        beg_k = 0
        end_k = nk

        if not (options.kpoint == None):
            print "# Warning: --kpoint option is discarded"
        
        if options.mode == None:
            beg_s = 0
            end_s = ns
        else:
            if len(options.mode.split(':')) == 1:
                beg_s = int(options.mode) - 1
                end_s = beg_s + 1
            elif len(options.mode.split(':')) == 2:
                arr = options.mode.split(':')
                beg_s, end_s = int(arr[0]) - 1, int(arr[1]) 
            else:
                sys.exit("Invalid usage of --mode for --calc=kappa_size")

        if options.length == None:
            max_len = 1000.0
            d_len = 1.0
        elif len(options.length.split(':')) == 2:
            arr = options.length.split(':')
            max_len, d_len = float(arr[0]), float(arr[1])
        else:
            sys.exit("Invalid usage of --length option")

    else:

        beg_k = 0
        end_k = nk

        if not (options.kpoint == None):
            print "# Warning: --kpoint option is discarded"

        if options.mode == None:
            beg_s = 0
            end_s = ns
        else:
            if len(options.mode.split(':')) == 1:
                beg_s = int(options.mode) - 1
                end_s = beg_s + 1
            elif len(options.mode.split(':')) == 2:
                arr = options.mode.split(':')
                beg_s, end_s = int(arr[0]) - 1, int(arr[1]) 
            else:
                sys.exit("Invalid usage of --mode for --calc=kappa_size")

    # Read Relaxation Times

    locate_tag("##Phonon Relaxation Time", fs)
    vel_tmp = [0]*3
    
    for i in range(nk):
        for j in range(ns):

            fs.readline()
            fs.readline()
            
            nk_equiv = int(fs.readline())
            
            # Read velocity

            for k in range(nk_equiv):

                data = fs.readline().split()
                vel_tmp[0], vel_tmp[1], vel_tmp[2] = [float(data[l]) for l in range(3)]

                for icrd in range(3):
                    for jcrd in range(3):
                        vel_sum[i, j, icrd, jcrd] += vel_tmp[icrd] * vel_tmp[jcrd]

                if k == 0:
                    for icrd in range(3):
                        vel_rep[i,j,icrd] = vel_tmp[icrd]

            for k in range(nt):
                tau[i, j, k] = float(fs.readline())
        
            fs.readline()
    
    # Calculate and write phonon properties

    if calc == "damping":
        print "# Temperature dependence of the damping function will be presented"
        print "# for phonon specified by kpoint", options.kpoint, " and mode", options.mode, "."
        print "# Temperature [k], Relaxation Time [ps], MFP [nm]"

        vel_norm = math.sqrt(vel_rep[target_k, target_s, 0]**2 + vel_rep[target_k, target_s, 1]**2 + vel_rep[target_k, target_s, 2]**2)
        
        for i in range(nt):
            print "%9.3f %15.7f %15.7f" % (temp[i], tau[target_k, target_s, i], tau[target_k, target_s, i]*vel_norm*0.001)

    elif calc == "tau":
        print "# Relaxation time at temperature ", temp[itemp], " K."
        print "# kpoint range ", beg_k + 1, end_k
        print "# mode   range ", beg_s + 1, end_s
        print "# ik, is, Frequency [cm^{-1}], Relaxation Time [ps], velocity [m/s],  MFP [nm]"

        for i in range(beg_k, end_k):
            for j in range(beg_s, end_s):
                vel_norm = math.sqrt(vel_rep[i,j,0]**2 + vel_rep[i,j,1]**2 + vel_rep[i,j,2]**2)
                print "%5i %5i %9.3f %15.7f %15.7f %15.7f" % (i + 1, j + 1, omega[i,j], tau[i,j,itemp], vel_norm, tau[i,j,itemp]*vel_norm*0.001)
                
    elif calc == "kappa":
        print "# Temperature dependence of thermal conductivity will be printed."
        print "# mode range", beg_s + 1, end_s
        print "# Temperature [K], kappa [W/mK] (xx, xy, xz, yx, yy, yz, zx, zy, zz)"

        kappa = np.zeros((nt, 3, 3))

        for itemp in range(nt):
            for i in range(nk):
                for j in range(beg_s, end_s):
                    Cv_tmp = Cv(omega[i,j],temp[itemp])
                    tau_tmp = tau[i,j,itemp]
                    for icrd in range(3):
                        for jcrd in range(3):
                            kappa[itemp,icrd,jcrd] += Cv_tmp*tau_tmp*vel_sum[i,j,icrd,jcrd]

        for itemp in range(nt):
            for icrd in range(3):
                for jcrd in range(3):
                    kappa[itemp,icrd,jcrd] *= factor

        for itemp in range(nt):
            print "%10.3F" % temp[itemp],
            for icrd in range(3):
                for jcrd in range(3):
                    print "%15.7g" % kappa[itemp, icrd, jcrd],
            print 
        
    elif calc == "kappa_size":
        print "# Size dependent thermal conductivity at temperature ", temp[itemp], " K."
        print "# L [nm], kappa [W/mK] (xx, xy, ...)"

        nlen = int(max_len / d_len)

        for ilen in range(nlen):
            kappa_size = np.zeros((3, 3))
            length = float(ilen) * d_len
            
            for k in range(nk):
                for j in range(beg_s, end_s):
                    tau_tmp = tau[k, j, itemp]
                    mfp_tmp = tau[k, j, itemp] * math.sqrt(vel_rep[k,j,0]**2 + vel_rep[k,j,1]**2 + vel_rep[k,j,2]**2) * 0.001
                    
                    if mfp_tmp <= length:
                        Cv_tmp = Cv(omega[i,j],temp[itemp])
                        for icrd in range(3):
                            for jcrd in range(3):
                                kappa_size[icrd,jcrd] += Cv_tmp*tau_tmp*vel_sum[k,j,icrd,jcrd]
                        
            for icrd in range(3):
                for jcrd in range(3):
                    kappa_size[icrd,jcrd] *= factor

            print "%15.7F" % length,
            for icrd in range(3):
                for jcrd in range(3):
                    print "%15.7g" % kappa_size[icrd, jcrd],
            print 
