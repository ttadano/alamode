#!/usr/bin/env python

from pylab import *
import sys
import matplotlib.pyplot as plt

il= 0
ip= 0
color= ['b', 'g', 'r', 'm', 'k', 'c', 'y', 'r']
lsty = ['-', '-', '-', '-', '--', '--', '--', '--']
psty = ['o']

for filename in sys.stdin:

    ifile= filename.rstrip('\n').split()

    if ifile[0] == 'l':
        f= open(ifile[1], 'r')
        k= []
        e= []
        e0= []

        for line in f:
            arrays = line.split()
            if arrays[0] != '#':
                k.append(float(arrays[0]))
                e0.append(float(arrays[1]))
                e.append(arrays[2:])

        plot(k,e0,linestyle=lsty[il],color=color[il],label=ifile[1])
        plot(k,e ,linestyle=lsty[il],color=color[il])
        il+= 1

    elif ifile[0] == 'p':
        f= open(ifile[1], 'r')
        k= []
        e= []
        e0= []

        for line in f:
            arrays = line.split()
            if arrays[0] != '#':
                k.append(float(arrays[0]))
                e0.append(float(arrays[1]))

        plot(k,e0,linestyle='None',
             marker=psty[0],color=color[ip],label=ifile[1])
        ip+= 1

    elif ifile[0] == 'kpath':
        kpind = []
        kpname= []
        kpind = ifile[1:]
        for i in range(len(kpind)):
            kpname.append(kpind[i])


    elif ifile[0] == 'kval':
        kpval= []
        for i in range(len(ifile[1:])):
            kpval.append(float(ifile[i+1]))

    else:
       break

ylabel("Frequency [cm${}^{-1}$]", fontsize=18)

xmin, xmax, ymin, ymax= plt.axis()

xmax = k[-1]
ymin = float(sys.argv[1])
ymax = float(sys.argv[2])
plt.axis([ xmin, xmax, ymin, ymax ])

plt.xticks(kpval, kpname, fontsize=18)
plt.yticks(fontsize=14)

plt.grid('on', which='major')

legend(loc='lower right')
plt.savefig('band_tmp.png', dpi=300, transparent=True)
show()


