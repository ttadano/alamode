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
* FFTW library
* `Eigen3 library <http://eigen.tuxfamily.org/>`_

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
VMD may be more useful to make an animation, but it may be replaced by any other visualization software which supports the XYZ format.

.. _XcrySDen : http://www.xcrysden.org
.. _VMD : http://www.ks.uiuc.edu/Research/vmd/

How to install
--------------

.. highlight:: bash

0. Install the LAPACK, MPI, FFTW, Boost C++, and Eigen3 libraries.

   To install the Boost C++ library, please download a source file from the `webpage <http://www.boost.org>`_ and
   unpack the file. Then, copy the 'boost' subdirectory to the include folder in the home directory (or anywhere you like).
   This can be done as follows::
    
    $ cd
    $ mkdir etc; cd etc
    (Download a source file and mv it to ~/etc)
    $ tar xvf boost_x_yy_z.tar.bz2
    $ cd ../
    $ mkdir include; cd include
    $ ln -s ../etc/boost_x_yy_z/boost .

  In this example, we made a symbolic link to the 'boost' subdirectory in ``$HOME/include``.
  Instead of installing from source, you can install the Boost library with `Homebrew <http://brew.sh>`_ on Mac OSX.

  In the same way, please install the Eigen3 include files as follows::

    $ cd
    $ mkdir etc; cd etc
    (Download a source file and mv it to ~/etc)
    $ tar xvf eigen-eigen-*.tar.bz2 (* is an array of letters and digits)
    $ cd ../
    $ cd include
    $ ln -s ../etc/eigen-eigen-*/Eigen .  

1. Download the package of ALAMODE from the download page or clone from the git repository.

2. Change directory to the location of the file and untar the file as follows::

	$ tar -xvzf alamode-x.y.z.tar.gz 

  This will create a directory alamode-x.y.z containing the following subdirectories:
  
  * alm/      : Source files for alm (force constant calculation)
  * anphon/   : Source files for anphon (phonon calculation)
  * external/ : Third-party include files
  * include/  : Commonly-used include files
  * tools/    : Small auxiliary programs and scripts
  * docs/     : Source files for making documents
  * example/  : Example files


.. highlight:: makefile


3. Edit the Makefiles

  In directories alm/ and anphon/, we provide sample Makefiles for Linux (with Intel compiler) and Mac OSX (with gcc). 
  Please copy one of them as ``Makefile`` and modify it appropriately for your environment.

  Here's a typical setting for Linux with Intel compiler::

    CXX = icpc 
    CXXFLAGS = -O2 -xHOST -openmp -std=c++11
    INCLUDE = -I../include -I$(HOME)/include

    CXXL = ${CXX}
    LDFLAGS = -mkl

    LAPACK = 
    LIBS = ${LAPACK}

  To enable OpenMP parallelization, please add the ``-openmp`` (Intel) or ``-fopenmp`` (gcc) option in ``CXXFLAGS``.
  In addition, the directory containing the boost/ and Eigen/ subdirectories must be given in ``INCLUDE``. 

4. Make executables by ``make`` command.

   If the compilation is successful, the binary file named **alm** (**anphon**) is created in the alm/ (anphon/) directory.
   To use some auxiliary scripts for post-processing and data conversion, please compile the programs in the tools directory as well.
   See README.md in the tools directory for details about the auxiliary programs.


