ALM: Output files 
-----------------

* ``PREFIX``.pattern_HARMONIC, ``PREFIX``.pattern_ANHARM?

 In these files, displacement patterns are printed in units of
 :math:`\boldsymbol{e}_{x,y,z}`.
 These files are created when ``MODE = suggest``.
 Patterns for anharmonic force constants are printed only when ``NORDER > 1``.

* ``PREFIX``.fcs

 Harmonic and anharmonic force constants in Rydberg atomic units.
 In the first section, only symmetry-reduced force constants are printed.
 All symmetry-related force constants are shown in the following section
 with the symmetry prefactor (:math:`\pm 1`).
 Created when ``MODE = fitting``.

* ``PREFIX``.xml

 An XML file containing the necessary information for performing
 phonon calculations.
 The files can be read by *anphon* using the ``FCSXML``-tag.
 Created when ``MODE = fitting``.


