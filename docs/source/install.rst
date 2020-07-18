Installation
============

Requirement
-----------

Mandatory requirements
~~~~~~~~~~~~~~~~~~~~~~

* C++ compiler (Intel compiler is recommended.)
* LAPACK library
* MPI library (OpenMPI, MPICH2, IntelMPI, etc.)
* `Boost C++ library <http://www.boost.org>`_
* FFTW3 library (not necessary when Intel MKL is available)
* `Eigen3 library <http://eigen.tuxfamily.org/>`_
* `spglib <https://atztogo.github.io/spglib/>`_

In addition to the above requirements, users have to get and install a first-principles package 
(such as VASP_, QUANTUM-ESPRESSO_, OpenMX_, or xTAPP_) or another force field package (such as
LAMMPS_) by themselves in order to compute harmonic and anharmonic force constants.

.. _VASP : http://www.vasp.at
.. _OpenMX : http://www.openmx-square.org
.. _QUANTUM-ESPRESSO : http://www.quantum-espresso.org
.. _xTAPP : http://frodo.wpi-aimr.tohoku.ac.jp/xtapp/index.html
.. _LAMMPS : http://lammps.sandia.gov


Optional requirements
~~~~~~~~~~~~~~~~~~~~~

* Python (> 2.6), Numpy, and Matplotlib
* XcrySDen_ or VMD_

We provide some small scripts written in Python for visualizing phonon dispersion relations, phonon DOSs, etc.
To use these scripts, one need to install the above Python packages.
Additionally, XcrySDen is necessary to visualize the normal mode directions and animate the normal mode.
VMD may be more useful to make an animation, but it may be replaced by any other visualization software which supports the XYZ format.

.. _XcrySDen : http://www.xcrysden.org
.. _VMD : http://www.ks.uiuc.edu/Research/vmd/

How to install
--------------

.. highlight:: bash

Here, we do not explain how to install a C++ compiler, LAPACK, MPI, and FFTW libraries because they are usually available on supercomputing systems.

Boost C++ and Eigen3 libraries (header files only)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some header files of Boost C++ and Eigen3 libraries are necessary to build ALAMODE binaries.
Here, we install header files of these libraries in ``$(HOME)/include``.
You can skip this part if these libraries are already installed on your system.
   
To install the Boost C++ library, please download a source file from the `webpage <http://www.boost.org>`_ and
unpack the file. Then, copy the 'boost' subdirectory to ``$(HOME)/include``.
This can be done as follows::
    
  $ cd
  $ mkdir etc; cd etc
  (Download a source file and mv it to ~/etc)
  $ tar xvf boost_x_yy_z.tar.bz2
  $ cd ../
  $ mkdir include; cd include
  $ ln -s ../etc/boost_x_yy_z/boost .

In this example, we place the boost files in ``$(HOME)/etc`` and create a symbolic link to the ``$(HOME)/boost_x_yy_z/boost`` in ``$(HOME)/include``.
Instead of installing from source, you can install the Boost library with `Homebrew <http://brew.sh>`_ on Mac OSX.

In the same way, please install the Eigen3 include files as follows::

  $ cd
  $ mkdir etc; cd etc
  (Download a source file and mv it to ~/etc)
  $ tar xvf eigen-eigen-*.tar.bz2 (* is an array of letters and digits)
  $ cd ../
  $ cd include
  $ ln -s ../etc/eigen-eigen-*/Eigen .  


If you have followed the instruction, you will see the following results::

  $ pwd
  /home/tadano/include
  $ ls -l
  total 0
  lrwxrwxrwx 1 tadano sim00 25 May 17  2017 boost -> ../etc/boost_1_64_0/boost
  lrwxrwxrwx 1 tadano sim00 38 May 17  2017 Eigen -> ../etc/eigen-eigen-67e894c6cd8f/Eigen/


spglib
~~~~~~

ALAMODE uses spglib to handle symmetries of crystal structures.
Please install it by following the instruction on the `spglib webpage <https://atztogo.github.io/spglib/install.html>`_.
Here, we assume spglib is installed in ``$(SPGLIB_ROOT)``.

Build ALAMODE
~~~~~~~~~~~~~

Download source
+++++++++++++++

From the download page::

  $ (visit https://sourceforge.net/projects/alamode/files/latest/download?source=files to download the latest version source)
  $ tar xvzf alamode-x.y.z.tar.gz
  $ cd alamode-x.y.z

From GitHub repository::

  $ git clone https://github.com/ttadano/alamode.git
  $ cd alamode
  $ git checkout master

The meaning of each subdirectory is as follows:

  * alm/      : Source files of alm (force constant calculator)
  * anphon/   : Source files of anphon (anharmonic phonon calculator)
  * external/ : Third-party include files
  * include/  : Commonly-used include files
  * tools/    : Small auxiliary programs and scripts
  * docs/     : Source files for making documents
  * example/  : Example files


Option 1. Build by Makefile
++++++++++++++++++++++++++++

ALAMODE contains two major codes, **alm** and **anphon**, and other small utility scripts.
In directories ``alm/``, ``anphon/``, and ``tools``, we provide sample Makefiles for Linux (Intel compiler) and Mac OSX (gcc, clang).
Please copy either of them, edit the options appropriately, and issue ``make`` command as follows::


    $ cd alm/
    $ cp Makefile.linux Makefile
    (Edit Makefile here)
    $ make -j

    $ cd ../anphon/
    $ cp Makefile.linux Makefile
    (Edit Makefile here)
    $ make -j

    $ cd ../tools/
    (Edit Makefile here)
    $ make -j

An example of the Makefiles is shown below:

.. literalinclude:: ../../alm/Makefile.linux
    :caption: **ALM Makefile.linux**
    :language: makefile
    :linenos:
    :lines: 7-17

The default options are expected to work with modern Intel compilers.

.. note::
   To build the binaries with the example Makefiles, you need to set ``SPGLIB_ROOT`` beforehand from the terminal as::

       $ export SPGLIB_ROOT=/path/to/spglib/installdir


Option 2. Build by CMake
++++++++++++++++++++++++

We also provide ``CMakeLists.txt`` for building binaries with CMake. 
To use this approach, you need  to install cmake version 3.1 or later.

Download source files from GitHub repository::

  $ git clone https://github.com/ttadano/alamode.git
  $ cd alamode
  $ git checkout develop

Then, build Makefiles with cmake::

  $ cd alamode
  $ mkdir _build; cd _build
  $ cmake -DSPGLIB_PATH=/path/to/spglib/installdir -DCMAKE_C_COMPILER=gcc-10 -DCMAKE_CXX_COMPILER=g++-10 ..

You can use ``-DCMAKE_C_COMPILER`` and ``-DCMAKE_CXX_COMPILER`` options to specify the compilers to build ALAMODE binaries. 
If these options are not given, cmake will detect the compilers automatically by referencing the environmental variables 
``${CC}`` and ``${CXX}``. 

The cmake options for popular compilers are shown below:

**Intel compiler**::

  $ cmake -DSPGLIB_PATH=/path/to/spglib/installdir -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_CXX_FLAGS="-O2 -xHOST" ..

**gcc**::

  $ cmake -DSPGLIB_PATH=/path/to/spglib/installdir -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_CXX_FLAGS="-O2 -march=native" ..

**gcc (Installed via Homebrew on macos)**::

  $ cmake -DSPGLIB_PATH=/path/to/spglib/installdir -DCMAKE_C_COMPILER=gcc-10 -DCMAKE_CXX_COMPILER=g++-10 -DCMAKE_CXX_FLAGS="-O2 -march=native" ..

You may need to change ``gcc-10`` (``g++-10``) to ``gcc-9`` (``g++-9``) or older versions installed on your system.

.. note::

    If cmake cannot find Boost or Eigen3 automatically, you need to give paths to the include directories by using ``-DBOOST_INCLUDE`` and ``-DEIGEN3_INCLUDE`` options as::
    
        $ cmake -DBOOST_INCLUDE=/path/to/boost/includedir -DEIGEN3_INCLUDE=/path/to/eigen3/includedir -DSPGLIB_PATH=/path/to/spglib/installdir -DCMAKE_C_COMPILER=gcc-10 -DCMAKE_CXX_COMPILER=g++-10  ..

Then, issuing
::

  $ make -j

will create all binaries in alm/, anphon/, and tools/ subdirectories under the current directory (_build). 
You can specify the binary to build, for example, as
::

  $ make alm -j

.. note::

    If the build of **anphon** fails due to an error related to FFTW3, you need to give the location of FFTW3 on your system by using the ``-DFFTW3_ROOT`` option, for example, as
    ::
      
      $ cmake -DFFTW3_ROOT=/path/to/fftw3/ -DSPGLIB_PATH=/path/to/spglib/installdir -DCMAKE_C_COMPILER=gcc-10 -DCMAKE_CXX_COMPILER=g++-10  ..

Set environmental variable ``LD_LIBRARY_PATH``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Finally, add the following line in ``$(HOME)/.bashrc`` (or ``.zshrc`` etc.)::

    (bash, zsh)
    export LD_LIBRARY_PATH=/path/to/spglib/installdir/lib:$LD_LIBRARY_PATH

    (csh, tcsh)
    setenv LD_LIBRARY_PATH /path/to/spglib/installdir/lib:$LD_LIBRARY_PATH

This is necessary when you link spglib dynamically.
