List of output files
--------------------

.. _reference_output:


Output files of alm
~~~~~~~~~~~~~~~~~~~

* ``PREFIX``.pattern_HARMONIC, ``PREFIX``.pattern_ANHARM?

 In these files, displacement patterns are printed in units of :math:`\boldsymbol{e}_{x,y,z}`.
 These files are created when ``MODE = suggest``. 
 Patterns for anharmonic force constants are printed only when ``NORDER > 1``.

* ``PREFIX``.fcs

 Harmonic and anharmonic force constants in Rydberg atomic units.
 In the first section, only symmetry-reduced force constants are printed.
 All symmetry-related force constants are shown in the following section with the symmetry prefactor (:math:`\pm 1`).
 Created when ``MODE = fitting``.

* ``PREFIX``.xml

 A XML file containing necessary information for performing phonon calculations.
 The files can be read by *anphon* using the ``FCSXML``-tag.
 Created when ``MODE = fitting``.


````

Output files of anphon
~~~~~~~~~~~~~~~~~~~~~~

.. |umulaut_u|    unicode:: U+00FC


* ``PREFIX``.bands

 Phonon dispersion along given :math:`k` paths in units of cm :sup:`-1`.
 Created when ``MODE = phonons`` with **KPMODE** = 1.

* ``PREFIX``.dos

 Phonon density of states (DOS). Atom projected phonon DOSs are also printed when ``PDOS = 1``.
 Created when ``MODE = phonons`` with **KPMODE** = 2.

* ``PREFIX``.tdos

 Two-phonon density of states for all irreducible :math:'k' points. 
 Created when ``MODE = phonons`` with **KPMODE** = 2 and ``TDOS = 1``.

* ``PREFIX``.thermo

 Constant volume heat capacity, vibrational entropy, internal energy, and vibrational free energy.
 Created when ``MODE = phonons`` with **KPMODE** = 2.

* ``PREFIX``.msd
 
 Mean-square-displacements of atoms.
 Created when ``MODE = phonons`` with **KPMODE** = 2 and ``PRINTMSD = 1``.

* ``PREFIX``.sps

 Total and mode-decomposed scattering phase space. 
 Created when ``MODE = phonons`` with **KPMODE** = 2 and ``SPS = 1``.

* ``PREFIX``.pr

 Participation ratio of every phonon modes. 
 Created when ``MODE = phonons`` and ``PRINTPR = 1``.

* ``PREFIX``.apr

 Atomic participation ratio of every phonon modes. 
 Created when ``MODE = phonons`` and ``PRINTPR = 1``.

* ``PREFIX``.phvel

 Phonon group velocity along given :math:`k` paths.
 Created when ``MODE = phonons`` with **KPMODE** = 1 and ``PRINTVEL = 1``.

* ``PREFIX``.phvel_all

 Magnitude of group velocity :math:`|\boldsymbol{v}|` of all phonon modes at the uniform :math:`k` grid. 
 Created when ``MODE = phonons`` with **KPMODE** = 2 and ``PRINTVEL = 1``.

* ``PREFIX``.evec

 Eigenvalues and eigenvectors of dynamical matrices.
 Eigenvalues are printed in Rydberg atomic units.
 Created when ``MODE = phonons`` with ``PRINTEVEC = 1``.

* ``PREFIX``.gru

 Gr\ |umulaut_u|\ neisen parameters along given :math:`k` paths.
 Created when ``MODE = phonons`` with **KPMODE** = 1 and ``GRUNEISEN = 1``.


* ``PREFIX``.gru_all

 Gr\ |umulaut_u|\ neisen parameters of all phonon modes at the uniform :math:`k` grid.
 Created when ``MODE = phonons`` with **KPMODE** = 2 and ``GRUNEISEN = 1``.

* ``PREFIX``.axsf

 Zone-center phonon modes with directions indicated by arrows.
 This file can be visualized by XcrySDen.
 Created when ``MODE = phonons`` with ``PRINTXSF = 1``.

* ``PREFIX``.anime???.axsf and ``PREFIX``.anime???.xyz

 Files for animating phonon modes. ??? is the mode number.
 Created when ``MODE = phonons`` with a proper ``ANIME``-tag.
 If ``ANIME_FORMAT = xsf``, axsf files will be created which can be displayed by XcrySDen.
 If ``ANIME_FORMAT = xyz``, xyz files will be created which can be visualized by VMD, Jmol, etc.

````

* ``PREFIX``.result

 In this file, phonon frequency, group velocity, and anharmonic phonon linewidths are printed.
 This file is updated during thermal conductivity calculations (``MODE = RTA``).
 In addition, this file is read when the restart mode is turned on (``RESTART = 1``).

* ``PREFIX``.kl

 Lattice thermal conductivity tensor.
 Created when ``MODE = RTA``.

* ``PREFIX``.kl_spec

 Spectra of lattice thermal conductivity. Only diagonal components will be saved.
 Created when ``MODE = RTA`` and ``KAPPA_SPEC = 1``.

* ``PREFIX``.gamma_isotope

 Phonon selfenergy due to isotope scatterings calculated by the Tamura's formula.
 Created when ``MODE = RTA`` and ``ISOTOPE = 2``.

````

* ``PREFIX``.scph_dymat

 Anharmonic dynamical matrix calculated on the :math:`k` grid defined by the ``KMESH_INTERPOLATE`` tag.
 This file is used to restart the SCPH calculation.

* ``PREFIX``.scph_bands

 Anharmonic phonon dispersion curves. The format is same as the ``PREFIX``.bands.

* ``PREFIX``.scph_fc2_correction

 This file contains :math:`\Delta D(\boldsymbol{q}) = D_{\mathrm{SCPH}}(\boldsymbol{q}) - D_{\mathrm{Harmonic}}(\boldsymbol{q})`.
 For the definition, see the :ref:`formalism of the SCPH calculation <formalism_SCPH>`.
