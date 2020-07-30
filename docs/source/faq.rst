.. _label_faq:

.. raw:: html

    <style> .red {color:red} </style>
    <style> .question {color:#cc0000; font-weight:bold} </style>


.. role:: red

.. role:: question


Frequently Asked Questions (FAQ)
================================

- :question:`The fitting error is very large (> 90%). Is it problematic? If so, how can I reduce the error?`
  
  The large fitting error can influence the accuracy of force constants. The most likely reason of the large fitting error is the non-zero residual forces in the original supercell structure (before making displacements of atoms). Even when the structure optimization is performed for a primitive cell with a relatively strict convergence criteria, the atomic forces in the supercell structure, which is generated from the primitive cell, may deviate from zero. To reduce the error associated with the residual forces, please use the ``--offset`` option of :red:`extract.py` when generating the displacement-force datasets. For example, in the case of VASP calculator, please issue
  ::

      $ python ${ALAMODE_ROOT}/tools/extract.py --VASP SPOSCAR --offset vasprun0.xml vasprun_harm*.xml > DFSET_harmonic

  Here, ``vasprun0.xml`` is the file obtained by running a VASP calculation for the original supercell (``SPOSCAR``).

  If the fitting error is still large after subtracting the offset components, please consider the following points:

  1. Please use ~15 decimal points for the fractional coordinates. For example, 1/3 should be 0.33333333333333 instead of 0.33333.

  2. Please check if the DFT calculation converges properly.

  3. Use a smaller displacement magnitude.


- :question:`How small the fitting error should be?`

  It depends on the Taylor expansion potential and displacement magnitude you choose. 
  
  In the standard harmonic calculation where ``--mag=0.01`` is used in :red:`displace.py` and the all harmonic interactions are considered (cutoff = None), the fitting error is usually less than 5% (~1--2% in most cases).

  In the calculation of third-order force constants with ``--mag=0.04``, the fitting error should be small as well. Indeed, in many cases, we obtain much smaller fitting errors (< 1%) than the harmonic case.

  In the temperature-dependent effective potential method, we try to fit the harmonic potential to the displacement-force datasets sampled by ab initio molecular dynamics at finite temperature. Therefore, the fitting error tends to be much larger (> 10%).

 
- :question:`What value should I use for the cutoff radius?`

  For the harmonic term, I would recommend using "None", which considers all harmonic interactions inside the supercell. This choice does not increase the computational costs that much because the number of displacement patterns does not change. Also, the harmonic dynamical matrix becomes exact at the commensurate q points only when the "None" option is selected. (Giving a very large cutoff radius has the same effect as giving "None".)

  For the anharmonic terms, you will need to increase the cutoff radii gradually and check the convergence of physical quantities, such as thermal conductivity and free energy, with respect to the cutoff values. In most cases, the cutoff radius of 10 Bohr is a good guess, but you may need to use a larger value for polar materials. So be careful.

- :question:`Why are the phonon dispersion curves discontinuous at the Brillouin zone boundaries?`

  Probably, you are wrongly using the supercell lattice vectors for the ``&cell`` field of the **anphon** code. If so, :red:`please use the primitive lattice vectors` for **anphon**. 

  .. Note::

      For the ``&cell`` field of **alm**, you need to give the supercell lattice vectors.


