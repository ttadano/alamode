Installation
============

Requirement
-----------

Mandatory requirements
~~~~~~~~~~~~~~~~~~~~~~

* C++ compiler (Intel compiler is recommended.)
* LAPACK library
* MPI library (Either OpenMPI, MPICH2, or IntelMPI)
* `Boost C++ library <http://www.boost.org>`_

In addition to the above requirements, users have to get and install a first-principles package 
(such as VASP_, Wien2k_, QUANTUM-ESPRESSO_, or xTAPP_) or another force field package (such as
LAMMPS_) by themselves in order to compute harmonic and anharmonic force constants.

.. _VASP : http://www.vasp.at
.. _Wien2k : http://www.wien2k.at
.. _QUANTUM-ESPRESSO : http://www.quantum-espresso.org
.. _xTAPP : http://frodo.wpi-aimr.tohoku.ac.jp/xtapp/index.html
.. _LAMMPS : http://lammps.sandia.gov


Optional requirements
~~~~~~~~~~~~~~~~~~~~~

* Python (> 2.6), Numpy, and Matplotlib
* XcrySDen_ or VMD_

We provide some small scripts written in Python (Python 2) for visualizing phonon dispersion relations, phonon DOSs, etc.
To use these scripts, one need to install the above Python packages.
Additionally, XcrySDen is necessary to visualize the normal mode directions and animate the normal mode.
VMD may be more useful to make an animation, but it may be replaced by any other visualization software which support the XYZ format.

.. _XcrySDen : http://www.xcrysden.org
.. _VMD : http://www.ks.uiuc.edu/Research/vmd/

How to install
--------------

.. highlight:: bash

1. Download the package from the download page or from the git repository.

2. Change directory to the location of the downloaded file and untar the file as follows::

	$ tar -xvzf alamode-x.y.z.tar.gz 

  This will create a directory alamode-x.y.z containing the following sub-directories:
  
  * alm/      : Source files for alm (force constant calculation)
  * anphon/   : Source files for anphon (phonon calculation)
  * external/ : Third-party include files
  * include/  : Commonly-used include files
  * tools/    : Small auxiliary programs and scripts
  * docs/     : Source files for making documents
  * example/  : Example files

3. Edit the Makefiles

  In directories alm/ and anphon/, we provide sample Makefiles for gcc and Intel compiler. 
  Please copy one of them as ``Makefile`` and modify it appropriately.
  To enable OpenMP parallelization, please add the ``-openmp`` (Intel) or ``-fopenmp`` (gcc) option in ``CXXFLAGS``.
  In addition, the directory containing the boost/ subdirectory must be given in ``INCLUDE``. 


4. Make executables by make command.

