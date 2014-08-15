Installation
============

Requirement
-----------

Mandatory requirements
~~~~~~~~~~~~~~~~~~~~~~

* C++ compiler (Intel compiler is highly recommended.)
* LAPACK library
* MPI library (Either OpenMPI, MPICH2, or IntelMPI)
* `Boost C++ library <http://www.boost.org>`_

In addition to the above requirements, users have to get and install a first-principles package 
(such as VASP_, Wien2k_, QUANTUM-ESPRESSO_, or xTAPP_) or another force field package (such as
LAMMPS_) by themselves in order to compute harmonic and anharmonic force constants.
Currently, the package provides some auxiliary tools for the VASP code, which 
may be useful to make POSCARs for displaced configurations and to extract atomic forces from
vasprun.xml files.

.. _VASP : http://www.vasp.at
.. _Wien2k : http://www.wien2k.at
.. _QUANTUM-ESPRESSO : http://www.quantum-espresso.org
.. _xTAPP : http://frodo.wpi-aimr.tohoku.ac.jp/xtapp/index.html
.. _LAMMPS : http://lammps.sandia.gov


Optional requirements
~~~~~~~~~~~~~~~~~~~~~

* Python, Numpy, and Matplotlib
* Xcrysden_

We provide some small scrips written in Python for visualizing phonon dispersion relations, phonon DOSs, etc.
To use these scripts, one need to install the above Python packages.
Xcrysden may be useful to visualize a zone-center normal mode. 

.. _Xcrysden : http://www.xcrysden.org

How to install
--------------

.. highlight:: bash

1. Download the package from the download page.

2. Change directory to the location of the downloaded file and untar the file as follows::

	$ tar -xvzf alamode-x.y.z.tar.gz 

  This will create a directory alamode-x.y.z containing the following sub-directories:
  
  * alm/ : Source files for alm (force constant calculation)
  * anphon/ : Source files for anphon (phonon calculation)
  * external/ : Third-party include files
  * include/ : Commonly-used include files
  * tools/ : Small auxiliary programs and scripts

3. Edit the Makefiles

  The directories alm/, anphon/, and tools/ contain separate Makefiles.
  Please modify the Makefiles appropriately by changing variables such as CXX, CXXFLAGS, INCLUDE.

4. Make executables by make command.