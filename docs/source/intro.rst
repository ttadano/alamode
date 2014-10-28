.. raw:: html

    <script type="text/javascript"><!--
      function gen_mail_to_link(name, lhs,rhs)
      {
       document.write("<br>");
       document.write(name + " (<A HREF=\"mailto");
       document.write(":" + lhs + "@");
       document.write(rhs + "\">" + lhs + "@" + rhs + "<\/A>)"); } 
       // --> </SCRIPT>


About
=====

What is ALAMODE ?
-----------------

**ALAMODE** stands for **A**\ nharmonic **LA**\ ttice **MODE**\ l, 
which is designed for estimating harmonic and anharmonic properties of lattice vibrations (phonons) in solids. 
The code is written in C++ and MPI/OpenMP parallelization is implemented.

Features
--------

ALAMODE consists of two main programs **alm** and **anphon** written in C++, and some subsidiary Python scripts.

The program *alm* can estimate harmonic and anharmonic interatomic force constants (IFCs) based on the supercell approach.
The program *anphon* can compute phonon properties using the estimated IFCs. Currently, it can compute the following quantities:

* Phonon dispersion, Phonon DOS, Atom-projected phonon DOS
* Thermodynamic functions (vibrational entropy, free energy, internal energy)
* Mean square displacement
* Heat capacity at constant volume
* Phonon linewidth due to 3-phonon interactions
* Phonon frequency shift
* Phonon self-energy due to phonon-isotope scatterings
* Lattice thermal conductivity

In addition, a python script may be used to estimate the cumulative thermal conductivity as a function of phonon mean-free-path.


Links
-----

* Download page  : http://sourceforge.net/projects/alamode
* Documentation  : http://alamode.readthedocs.org
* Git repository : http://github.com/ttadano/alamode


License
-------

.. |copy|   unicode:: U+000A9 

Copyright |copy| 2014 Terumasa Tadano. See the LICENSE.txt file for license
rights and limitations (MIT). 


If you used ALAMODE, please cite the following article:

  T\. Tadano, Y. Gohda, and S. Tsuneyuki, J. Phys.: Condens. Matter **26**\ , 225402 (2014) [Link_].

.. _Link : http://iopscience.iop.org/0953-8984/26/22/225402/


Acknowledgement
---------------

This project is support by a Grant-in-Aid for Scientific Research on Innovative Areas 
'Materials Design through Computics: Complex Correlation and Non-Equilibrium Dynamics'.
(http://computics-material.jp)


Author
------


.. raw:: html

    <script>gen_mail_to_link('Terumasa TADANO', 'terumasa.tadano','gmail.com')</script>

Current affiliation: Department of Physics, The University of Tokyo, Japan

