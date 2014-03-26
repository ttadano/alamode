from scipy import *
import numpy as np
from scipy import integrate
import sys

kelvin_to_Ryd = 6.33363081e-6
cminv_to_Ryd = 9.11267002e-6

dosfile = sys.argv[1]

temp_start = float(sys.argv[2])
temp_end = float(sys.argv[3])
temp_step = float(sys.argv[4])

temp = np.arange(temp_start, temp_end, temp_step)
tmp = loadtxt(dosfile)

omega = tmp[:,0]*cminv_to_Ryd
for i in range(len(omega)):
    if omega[i] < 1.0e-12:
        omega[i] = 1.0e-12

dos = tmp[:,1]

F = zeros(len(temp))

for i in range(len(temp)):
    x = temp[i]
    x *= kelvin_to_Ryd

    if x < 1.0e-10:
        x = 1.0e-10

    F_local = dos*(0.5*omega + x*log(1.0-exp(-omega/x)))
    F[i] = integrate.simps(F_local, omega)

    print "%10.5f %15.7e" % (temp[i], F[i])
