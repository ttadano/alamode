# Ver. 0.9.3 (2014-12-12)

## New

- New tag FC3XML (alm)

- New tag SPS (anphon) for calculating the scattering phase space  

- New mode calc=kappa_boundary in analyze_phonons.py

## Changes

- Changed the filename of TDOS from \*.dos2 to \*.tdos

## Fix

- Fixed a bug in phonon DOS and two-phonon DOS

- Fixed a bug in anphon regarding the NSYM variable

- Fixed a bug related to the isotope scattering rate

- Fixed xml_writer to be compatible with boost 1.56.0

# Ver. 0.9.2 (2014-10-29)

## New

- Python auxiliary tool for Quantum-ESPRESSO

- Tutorial added

- Added displace.py and extract.py in tools/ directory

## Changes

- Different algorithm for anharmonic scattering rates (a little performance improvement is expected)

- PRINTVEL = 1 now print the xyz component

- Renamed phonon_thermodynamics.{cpp,h} as thermodynamics.{cpp,h}

## Fix

- Fixed a bug in analyze_phonons.cpp

- Fixed a MPI-related bug in relaxationc.pp

- Avoid unnecessary memory allocation in relaxation.cpp

- Fixed a bug in the rank-revealing algorithm of alm

# Ver. 0.9.1 (2014-08-24)

## New

- Python auxiliary tool for xTAPP

- alm now checks if DFILE and FFILE have enough lines

## Changes

- Change the format of eigenvector *.evec file

- Tmin and Tmax are now inclusive


## Fix

- Fixed many minor bugs

# Ver. 0.9.0 (2014-08-15)

-  First release of ALAMODE
