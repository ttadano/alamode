#!/usr/bin/env python
"""

Script for generating animation video of
temperature change obtained from MD simulations.

Usage:animete_T.py filename npoint ymin, ymax

"""
from scipy import *
import matplotlib.pyplot as plt
import numpy as np
import sys, os

if __name__ == '__main__':

    file_T = sys.argv[1]
    npoint = int(sys.argv[2])
    ymin = float(sys.argv[3])
    ymax = float(sys.argv[4])

    T = loadtxt(file_T)

    xmin = 0
    xmax = T[npoint-1,0]

    nrow, ncol = np.shape(T)

    ndata = nrow / npoint
    files = []

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for i in range(ndata):
        ax.cla()
        ax.axis([xmin, xmax, ymin, ymax])
        ax.plot(T[i*npoint:(i+1)*npoint - 1,0], T[i*npoint:(i+1)*npoint - 1,1])
        fname = '_tmp%04d.png' % i
        print 'Saving frame', fname
        fig.savefig(fname)
        files.append(fname)

    print 'Making movie animation.mpg - this make take a while'
    os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=10 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o animation.mpg")

    print
    print 'Removing temporary png files'

    for i in files:
        command = "rm -rf " + i
        os.system(command)

    print 'Done'
