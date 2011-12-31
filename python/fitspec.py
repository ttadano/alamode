#!/usr/bin/env python
from scipy import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import sys

respower = lambda v, x, y: (powerfunc(v,x) - y)
respower2 = lambda v, x, y: ((powerfunc(v,x) - y)*sqrt(y))

def powerfunc(v,x):
    func = 0.0
    nm = (len(v) - 1) / 2
    for i in range(nm):
        j = 2*i + 1
        vtmp = v[j+1]
#        func +=  vtmp*( 1.0 / (1.0 + (vtmp*(x + v[j]))**2) + 1.0 / (1.0 + (vtmp*(x - v[j]))**2))
        func +=  vtmp / (1.0 + (vtmp*(x - v[j]))**2)
    return func*v[0]

infile = sys.argv[1]
powerspec = loadtxt(infile)

cinv_to_ps = 0.01 / 299792458 * 1.0e+12

try:
    nrow, ncol = shape(powerspec)
except:
    nrow = size(signal)
    ncol = 1
    powerspec = np.reshape(powerspec, (nrow,ncol))

print "#", nrow, ncol

datafile = open(infile, 'r')
tmp = datafile.readline().split()
nmode = int(tmp[1])
print "# nmode =", nmode

nspec = ncol - 2
xpeaks = zeros((nspec, nmode))

xmax = 0.0

# peaks for each spectrum
nmode_k = []

for i in range(nspec):
    line = datafile.readline().split()
    line = datafile.readline().split()
    nmode_k.append(int(line[1]))

    for j in range(nmode_k[i]):
        line = datafile.readline().split()
        xpeaks[i,j] = float(line[1])

    xmax = max(max(xpeaks[i,:]), xmax)

xmax = xmax*1.1
x = powerspec[:,0]

print "# k-points"
print "# Frequency[cm^-1] with err.     Life Times [ps] with err. (2*nmode data)"
print "#"

xfitted = zeros((nrow, nspec))

for i in range(nspec):

    params = []
    params.append(0.1)
    for j in range(nmode_k[i]):
        params.append(xpeaks[i,j])
        params.append(1.0)

    y = powerspec[:,i+2]
    param_out, cov_out, info, mesg, success = leastsq(respower2, params, args=(x,y), full_output=1)

    chisq = (info["fvec"]**2).sum()
    dof = len(y) - len(params)
    if cov_out != None:
        cov_out = cov_out * chisq / dof
    
    print "# k-point %7i" % (i+1)
    for j in range(nmode_k[i]):
        k = 2*j + 1
        if cov_out != None:
            sys.stdout.write("# %8.5f %8.5f   " % (param_out[k], sqrt(cov_out[k][k])))
            sys.stdout.write("%8.5f %8.5f" % (param_out[k+1], sqrt(cov_out[k+1][k+1])))
        else:
            sys.stdout.write("# %8.5f, None" % param_out[k])
            sys.stdout.write("%8.5f, None" % param_out[k+1])

        sys.stdout.write("\n")
        
    if int(sys.argv[2]) == 1:
        plt.plot(x,y)
        plt.plot(x,powerfunc(param_out,x))
    else:
        xfitted[:,i] = powerfunc(param_out, x[:])

if int(sys.argv[2]) == 1:
    plt.yscale('log')
#    plt.axis([-xmax, xmax, 1.0e-6, 1])
    plt.axis([0, xmax, 1.0e-8, 1])
    plt.show()
else:
    for j in range(nrow):
        sys.stdout.write("%15.7e" % x[j])
        for i in range(nspec):
            sys.stdout.write("%15.7e" % xfitted[j,i])
            sys.stdout.write("%15.7e" % (powerspec[j,i+2] - xfitted[j,i]))
        sys.stdout.write('\n')
