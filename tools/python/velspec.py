#!/usr/bin/env python

from scipy import *
from spectrum import *
import numpy as np
import sys

infile = sys.argv[1]
timestep = float(sys.argv[2])
signal = loadtxt(infile)
fs = 1.0e-15

datafile = open(infile, 'r')

freqlist = []

for i in range(4):
    line = datafile.readline()

freqlist = datafile.readline().split()
freqlist.remove('#')

try:
    nrow, ncol = shape(signal)
except:
    print "Input should be like (real.v[1] complex.v[1] ... real.v[N] complex.v[N])"
    exit

ncol = ncol / 2
nmode = 1

try:
    if sys.argv[3] == "acorr":
        usecorr = True
    elif sys.argv[3] == "direct":
        usecorr = False
    else:
        print "3rd input should be \"acorr\" or \"direct\""
        exit
except:
    print "3rd parameter is not given. \"direct\" method will be used."
    usecorr = False

try:
    if sys.argv[4] == "sum":
        sumflag = True
        nmode = int(sys.argv[5])
        print "#", nmode

        for i in range(ncol/nmode):
            print "# k-point", i+1
            xpeak = []
            xpeak.append(float(freqlist[nmode*i]))
            for j in range(nmode-1):
                k = j + 1
                tmp = float(freqlist[nmode*i+k])
                if tmp != xpeak[len(xpeak)-1]:
                    xpeak.append(tmp)
            print "# %i" % len(xpeak)
            for j in range(len(xpeak)):
                print "# %15.7f" % xpeak[j]

    elif sys.argv[4] == "raw":
        sumflag = False
    else:
        print "4th input should be \"sum\" or \"raw\""
        exit
except:
    print "#4th parameter is not given. \"sum\" will be used."
    sumflag = True


try:
    if sys.argv[6] == "minus":
        print_minus = True
    elif sys.argv[6] == "plus":
        print_minus = False
    else:
        print "4th input should be \"plus\" or \"minus\""
        exit
except:
    print_minus = False

sinv_to_cinv = 0.01 / 299792458.0
sinv_to_THz = 1.0e-12

nk = ncol/nmode

#print "#", nrow, ncol
#print "# sampling frequency", sinv_to_cinv / (fs*timestep)

sigcomplex = zeros((nrow, ncol), dtype=complex128)

for i in range(ncol):
    sigcomplex[:,i] = signal[:,2*i] + signal[:,2*i+1]*1j

nf = nrow
freq = np.fft.fftfreq(nrow, d = timestep*fs)
freq2 = freq * sinv_to_THz
freq = freq * sinv_to_cinv

powerspec = zeros((nf, ncol))

for i in range(ncol):
    avg = np.average(sigcomplex[:,i])
    if usecorr:
        sigcomplex[:,i] = sigcomplex[:,i] - avg

for i in range(ncol):
    if usecorr:
        ar, p, k = aryule(sigcomplex[:,i], 100)
#        ar, rho, ref = arburg(sigcomplex[:,i], 100)
        powerspec[:,i] = arma2psd(ar, NPSD=nrow)
#        fourier = np.correlate(sigcomplex[:,i],conjugate(sigcomplex[:,i]), mode='full') / float(nf)
#        powerspec[:,i] = np.fft.fft(fourier,n=nrow)
    else:
        sigcomplex[:,i] *= hamming(nrow)
        fourier = np.fft.fft(sigcomplex[:,i],n=nrow)
        powerspec[:,i] = (fourier[:]*conjugate(fourier[:])).real

    powerspec[:,i] = powerspec[:,i] / sum(powerspec[:,i])

index_sort = np.argsort(freq, kind = 'quicksort')

if sumflag:
    powerspec_sum = zeros((nf, nk))
    for i in range(nk):
        for j in range(nmode):
            powerspec_sum[:,i] += powerspec[:,i*nmode+j]
    for i in range(nf-1):
        if (not print_minus) and (freq[index_sort[i]] < 0.0) : continue
        sys.stdout.write("%15.7e %15.7e" % (freq[index_sort[i]], freq2[index_sort[i]]))
        for j in range(nk):
            sys.stdout.write("%15.7e" % powerspec_sum[index_sort[i],j])
        sys.stdout.write("\n")
else:
    for i in range(nf-1):
        sys.stdout.write("%15.7e %15.7e" % (freq[index_sort[i]], freq2[index_sort[i]]))
        for j in range(ncol):
            sys.stdout.write("%15.7e" % powerspec[index_sort[i],j])
        sys.stdout.write("\n")
