#!/usr/bin/env python

from scipy import *
import sys
import matplotlib.pyplot as plt

file = sys.argv[1]
data = loadtxt(file)

plt.fill(data[:,0], data[:,1], 'b', alpha = 0.5)
plt.xlabel("Energy [cm${}^{-1}$]", fontsize=18)
plt.ylabel("Phonon DOS", fontsize=18)
plt.grid("on")
plt.show()


