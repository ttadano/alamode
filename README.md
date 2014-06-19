# ALAMODE
* Version 0.9.0 (Beta)

- - -

## Introduction
Program package ALAMODE is designed for estimating anharmonic force constants of solids 
based on the supercell approach and subsequent calculations of anharmonic phonon properties, 
such as Gruneisen parameter, phonon self-energy and lattice thermal conductivity.

The package mainly consists of the following 2 programs:

* **alm**: Extract harmonic and anharmonic force constants from first-principles calculations
* **anphon** : Compute anharmonic phonon properties

## Prerequisite
* C++ compiler
* LAPACK libarary
* MPI library
* Boost C++ library

## Install
The directories alm/, anphon/, and tools/ contain separate Makefiles.
Please modify the Makefiles appropriately by changing variables such as 
CXX, CXXFLAGS, or MPICXX. Then, execute "make" will create the binary for
each program.
More detailed instruction may be found in the user's guide.


## Documentation
Please refer to the following webpage.
http://alamode.readthedocs.org


## License
Copyright (c) 2014 Terumasa Tadano
This software is released under the MIT license. 
For license rights and limitations, see LICENSE.txt file.

## Author
Terumasa Tadano (terumasa.tadano{at}.gmail.com)
