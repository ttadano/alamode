# Ver. 1.0.0 (2017-11-21)

## New

- Self-consistent phonon calculation (``MODE = SCPH``) 
- Ewald summation for the non-analytic term of the dynamical matrix (``NONANALYTIC = 3``).
- Support of the ``CLASSICAL`` option.
- Python auxiliary script for LAMMPS

## Changes

- P+ and P- are printed seperately in PREFIX.sps when ``SPS = 1``
- Use C++11 standard. From this version, the C++ compiler must support the C++11 standard.
- The **anphon** code symmetrize the Born effective charges

## Fix

- Loosen the tolerance to detect the multiplicity of force constants. This is an important fix for low symmetry structures.
- Fixed a problem of the restart mode of ``MODE = RTA``


# Ver. 0.9.8 (2016-7-14)

## New

- New tag ``HESSIAN`` in program **alm** for printing entire Hessian matrix
- New tag ``KAPPA_SPEC`` in **anphon** for calculating spectra of thermal conductivity
- New option ``--offset`` in **extract.py** for subtracting residual forces (displacements) in an equilibrium structure from training data sets. We recommended to use this option if internal coordinates of atoms (Wyckoff potisions) have free parameters.
- New option ``--isotope`` in **analyze_phonons.py**

## Changes

- Improve the performance of thermal conductivity calculations with the tetrahedron method (``ISMEAR=-1``). The new version is <Font color='red'>more than 3 times faster</font> than the previous version.
- Improve the efficiency of the algorithm for generating constraings for the translational invariance
- Avoid 'NaN' in thermal conductivity calculations with imaginary branches
- Loosen the default value of ``TOLERANCE`` for detecting crystal symmetry
- Stop printing the CLASSICAL entry in PREFIX.result files
- Change the unit of the smearing width written in PREFIX.result files (a.u. --> cm^-1)

## Fix

- Fixed an issue regarding the phonon-isotope scattering rate
- Fixed a bug in relaxation.cpp regarding the permutation symmetry of q mesh
- Fixed an invalid memory reference in **analyze_phonons.cpp**
- Added a routine to check the consistency of the crystal structure when the ``FC2XML`` tag is used in **anphon**


# Ver. 0.9.7 (2016-3-10)

## New

- New option ``SPS = 2`` in anphon
- New tag ``MAGMOM`` for considering collinear spin (alm). The format is same as VASP.
- ALAMODE logo 

## Changes

- From this version, the variable 'NNP' for the translational part of 
  space group operations will no longer be used. Because of this, 
  file formats of SYMM_INFO and SYMM_INFO_PRIM also changed.

- Comment lines in the &cell, &position, and &kpoints fields are supported.

- Any words starting from 'N' or 'n' can be used instead of 'None' for the &cutoff field.

## Fix

- Fixed a bug in fitting.cpp (alm) regarding out-of-range of an array.
- Fixed a bug in the calculation of the maximum distance for multi-body interaction clusters.
- Fixed a bug in the calculation of constraints for translational invariance. 
  This bug appeared in previous versions with some specific structures/cutoff radii.
- anphon code incorrectly printed RMSD instead of MSD in previous versions. Fixed.
- qe2alm.cpp now works with ibrav=0.
- Fixed several issues in displace.py and extract.py for processing QE files.
- ``PRINTPR`` now works properly when KPMODE = 0.


# Ver. 0.9.6 (2015-09-14)

## New

- New option ``ICONST=11`` in **alm** which considers translational invariance algebraically

- ``FC2XML``-tag in **anphon**. Using this tag, it is now possible to employ different size of 
   supercells for harmonic and anharmonic terms.

- qe2alm.cpp in the tools/ directory which converts Quantum-ESPRESSO fc files to ALAMODE xml files.

## Changes

- Changed the meaning of ``--calc=cumulative`` in analyze_phonons.py. 
   The ``--calc=cumulative`` in the previous version is now equivalent to ``--calc=cumulative2``.

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

- Implement ``ISOTOPE = 2`` to print the selfenergy due to phonon-isotope scatterings (anphon)

## Changes

- Stop printing "nzero" in alm because it may be confusing.

## Fix

- Fixed a bug related to PDOS (anphon)

- Fixed issue in displace.py for QE


# Ver. 0.9.3 (2014-12-12)

## New

- New tag ``FC3XML`` (alm)

- New tag ``SPS`` (anphon) for calculating the scattering phase space  

- New mode ``--calc=kappa_boundary`` in analyze_phonons.py

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

- ``PRINTVEL = 1`` now print the xyz component

- Renamed phonon_thermodynamics.{cpp,h} as thermodynamics.{cpp,h}

## Fix

- Fixed a bug in analyze_phonons.cpp

- Fixed a MPI-related bug in relaxationc.pp

- Avoid unnecessary memory allocation in relaxation.cpp

- Fixed a bug in the rank-revealing algorithm of alm

# Ver. 0.9.1 (2014-08-24)

## New

- Python auxiliary tool for xTAPP

- alm now checks if ``DFILE`` and ``FFILE`` have enough lines

## Changes

- Change the format of eigenvector *.evec file

- Tmin and Tmax are now inclusive


## Fix

- Fixed many minor bugs

# Ver. 0.9.0 (2014-08-15)

-  First release of ALAMODE
