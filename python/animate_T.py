#!/usr/bin/env python
"""

Script for generating animation video of
temperature change obtained from MD simulations.

Usage:animete_T.py filename npoint ymin ymax fps

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
    fps = sys.argv[5]
    dt = float(sys.argv[4])

    ifs = open(file_T)
    ipoint = 0
    iframe = 0
    x = []
    temperature = []
    files = []

    fig = plt.figure()
    ax = fig.add_subplot(111)

    xmin = 0
    xmax = None

    for line in ifs:
        if line.strip().startswith("#"):
            continue
        
        ipoint += 1
        a, b = line.split()
        x.append(float(a))
        temperature.append(float(b))

        if ipoint == npoint:
            if xmax == None:
                xmax = x[npoint - 1]

            ax.cla()
            ax.axis([xmin, xmax, ymin, ymax])
            ax.plot(x, temperature)
            time_str = 't = ' + str((iframe + 1) * dt) + ' ps'
            plt.text((xmax - xmin)*0.75, (ymax - ymin)*0.8, time_str)
            fname = '_tmp%05d.png' % iframe
            print 'Saving frame', fname
            fig.savefig(fname)
            files.append(fname)
            iframe += 1
            ipoint = 0

            # Clear history
            x = []
            temperature = []

    print 'Making movie animation.mpg - this make take a while'
    command = "mencoder 'mf://_tmp*.png' -mf type=png:fps=" + fps + " -ovc lavc -lavcopts vcodec=mpeg4:threads=4:vbitrate=1500 -oac copy -o animation.mp4"
    os.system(command)

    print
    print 'Removing temporary png files'

    for i in files:
        command = "rm -rf " + i
        os.system(command)

    print 'Done'
