#!/usr/bin/env python

from math import sqrt

def val_and_error(x, y, sigx, sigy):
    return x/y, sqrt(sigx**2+(x/y)**2*sigy**2)/y

message1 = 'Heat Flux J (and its error) in unit of nW/(AA^2):\n'
message2 = 'Temperature Gradient (and its error) in unit of K/nm:\n'

J, error_J  = map(float, raw_input(message1).split())
dT, error_dT= map(float, raw_input(message2).split())

J *= 100.0
error_J *= 100.0

print 'Thermal Conductivity (and error) in unit of W/mK\n %10f %10f' \
      % (val_and_error(J, dT, error_J, error_dT))

print 'Inverse Thermal Conductivity (and error)\n %10e %10e' \
      % (val_and_error(dT, J, error_dT, error_J))

