#!/usr/bin/env python

from pylab import *
import sys
import matplotlib.pyplot as plt
import numpy as np

energy = []
dos = []

file = sys.argv[1]
data = loadtxt(file)

plt.fill(data[:,0], data[:,1], 'b', alpha=0.5)
plt.xlabel("Energy [cm${}^{-1}$]", fontsize=18)
plt.ylabel("Phonon DOS", fontsize=18)
show()
