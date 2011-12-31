#!/usr/bin/env python
from scipy import *
import numpy as np
import sys

infile = sys.argv[1]
npar = int(sys.argv[2])
prefix = sys.argv[3]
objfile = open(infile)
line = objfile.readline().split()
line.remove("#")

nk, nmode, nrow, dt = line
nk = int(nk)
nmode = int(nmode)
nrow = int(nrow)
dt = float(dt)

objfile.readline()
objfile.readline()
objfile.readline()

line = objfile.readline().split()
line.remove("#")
omega = array(line, dtype=float64)

nrow_new = int(2**round(log(float(nrow))/log(2.0) + 0.5))
nf = nrow_new / 2 + 1
freq = np.fft.fftfreq(nrow_new, d = dt)

if (nk*nmode)%npar != 0 or nk*nmode < npar:
    print "invalid npar"
    exit

nset = nk*nmode / npar

powerspec = zeros((nf, npar))
acorr = zeros((nrow_new, npar))

for i in range(nset):
    emode = loadtxt(infile, usecols=range(i*npar, (i+1)*npar))

    if npar == 1:
        emode = emode.reshape(len(emode),1)
    
    for j in range(npar):
        avg = np.average(emode[:,j])
        emode[:,j] = emode[:,j] - avg
    
    for j in range(npar):
        fourier = np.fft.rfft(emode[:,j], n = nrow_new)
        powerspec[:,j] = fourier[:]*conjugate(fourier[:])
        acorr[:,j] = np.fft.irfft(powerspec[:,j])
        acorr[:,j] = 2.0*acorr[:,j]/float(nrow_new)
        acorr[:,j] = acorr[:,j] / acorr[0,j]
    
    outname = prefix + ".acorr"
    
    ofile = open(outname, 'w')
    ofile.write("#             ")
    for j in range(i*npar, (i+1)*npar):
        ofile.write("%15.7f" % omega[j])
    ofile.write("\n")
    for j in range(nrow/2):
        ofile.write("%15.7e" % (float(j)*dt*1.0e+12))
        for k in range(npar):
            ofile.write("%15.7e" % acorr[j,k])
        ofile.write("\n")

    ofile.close()


