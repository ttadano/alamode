#!/usr/bin/env python

import sys
from numpy import *
import numpy as np

file = open(sys.argv[1], 'r')
line = file.readline()

nat, natcell, nstep = map(int, line.split())

nstart = natcell*(int(sys.argv[2]) - 1)
nend   = natcell*int(sys.argv[3])
nat2   = nend - nstart
nstep  = 200

tmp = []

print '# Start reading velocity data'

for t in range(nstep):
    for i in range(nat):
        line = file.readline()
        if nstart <= i < nend:
            tmp.append(line.split())

vel_tmp = np.array(tmp, dtype=float)
vel     = np.reshape(vel_tmp,(nstep,nat2,3))

del vel_tmp, tmp

Z    = np.zeros((nstep), float)

print '# Start generating autocorrelation function'

for m in range(nstep):
    for n in range(nstep - m):
        for i in range(nat2):
            Z[m] += vel[n+m][i][2]*vel[n][i][2]
    Z[m] /= float(nstep-m)

print '# Start FFT'

freq = np.fft.rfft(Z)

for m in range(nstep/2+1):
    print "%d %-10e" % (m, np.abs(freq[m]))

file.close()
