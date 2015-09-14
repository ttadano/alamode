# Ver. 0.9.6 (2015-09-14)

## New

- New option ICONST=11 for ALM which considers translational invariance algebraically

- FC2XML-tag for ANPHON. Using this tag, it is now possible to employ different size of 
   supercells for harmonic and anharmonic terms.

- qe2alm.cpp in the tools/ directory which converts Quantum-ESPRESSO fc files to ALAMODE xml files.

## Changes

- Changed the meaning of --calc=cumulative in analyze_phonons.py. 
   The --calc=cumulative in the previous version is now equivalent to --calc=cumulative2.

- Updated tutorial

## Fix

- Fix bugs related to memory allocation

- Fix an issue in extract.py for QE


# Ver. 0.9.5 (2015-06-28)

## New

- PERIODIC tag for low dimensional systems

## Changes

- Print distances for anharmonic terms in PREFIX.fcs

- Improved the way to define the multiplicity of force constants

- OpenMP parallelization for generating constraints between force constants

## Fix

- Fixed a minor bug in interaction.cpp

- Fixed a bug in displace.py for QE

# Ver. 0.9.4 (2015-02-16)

## New

- New tag PRINTPR (anphon) for calculating (atomic) participation ratio

- Implement ISOTOPE = 2 to print the selfenergy due to phonon-isotope scatterings (anphon)

## Changes

- Stop printing "nzero" in alm because it may be confusing.

## Fix

- Fixed a bug related to PDOS (anphon)

- Fixed issue in displace.py for QE


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
