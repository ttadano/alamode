
Making input files for *alm*
----------------------------

.. _reference_input_alm:

Format of input files
~~~~~~~~~~~~~~~~~~~~~

Each input file should consist of entry fields.
Available entry fields are 

**&general**, **&interaction**, **&cutoff**, **&cell**, **&position**, and **&fitting**.


Each entry field starts from the key label **&field** and ends at the terminate character "/". (This is equivalent to Fortran namelist.) 

For example, &general entry field of program *alm* should be given like

::

  &general
    # Comment line
    PREFIX = prefix
    MODE = fitting
  /

Multiple entries can be put in a single line. Also, characters put on the right of sharp (“#”) will be neglected. Therefore, the above example is equivalent to the following::
  
  &general
    PREFIX = prefix; MODE = fitting  # Comment line
  /

Each variable should be written inside the appropriate entry field.


.. _label_inputvar_alm:

List of input variables
~~~~~~~~~~~~~~~~~~~~~~~


"&general"-field
++++++++++++++++


* **PREFIX**-tag : Job prefix to be used for names of output files

 :Default:  None
 :Type: String

````

* **MODE**-tag = fitting | suggest

 ========= ===========================================================
  fitting  | Perform fittings to estimate harmonic and anharmonic IFCs. 
           | This mode requires appropriate DFILE and FFILE.

  suggest  | This mode suggests the displacement patterns necessary 
           | to estimate harmonic and anharmonic IFCS.
 ========= ===========================================================

 :Default: None
 :Type: String

````

* **NAT**-tag : Number of atoms in the supercell

 :Default: None
 :Type: Integer

````

* **NKD**-tag : Number of atomic species

 :Default: None
 :Type: Integer

````

* **KD**-tag = Name[1], ... , Name[``NKD``]

 :Default: None
 :Type: Array of strings
 :Example: In the case of GaAs with ``NKD = 2``, it should be ``KD = Ga As``.

````

* NSYM-tag = 0 | 1 | nsym

 ===== ==========================================================
   0   | The program automatically generates the crystal symmetry 
       | operations (rotational and translational parts). 
       | When ``PRINTSYM = 1``, symmetry operations will be saved 
       | in the file “SYMM_INFO”.

   1   | Only the identity operation will be considered.
  nsym | "nsym" symmetry operations will be read from "SYMM_INFO" 
       | file.
 ===== ==========================================================

 :Default: 0
 :Type: Integer

````

* TOLERANCE-tag : Tolerance for finding symmetry operations
  
 :Default: 1.0e-6
 :Type: Double

````

* PRINTSYM-tag = 0 | 1

 === ====================================================
  0   Symmetry operations won’t be saved in “SYMM_INFO”
  1   Symmetry operations will be saved in “SYMM_INFO”
 === ====================================================

 :Default: 0
 :type: Integer

````

* PERIODIC-tag = PERIODIC[1], PERIODIC[2], PERIODIC[3] 

 ===== ====================================================
   0   | Do not consider periodic boundary conditions when
       | searching for interacting atoms.

   1   | Consider periodic boundary conditions when
       | searching for interacting atoms.
 ===== ====================================================

 :Default: 1 1 1
 :type: Array of integers
 :Description: This tag is useful for generating interacting atoms in low dimensional systems. When ``PERIODIC[i]`` is zero, periodic boundary condition is turned off along the direction of the lattice vector :math:`\boldsymbol{a}_{i}`.

````

* HESSIAN-tag = 0 | 1

 ===== =====================================================================
   0    Do not save the Hessian matrix
   1    Save the entire Hessian matrix of the supercell as PREFIX.hessian.
 ===== =====================================================================

 :Default: 0
 :type: Integer

````

"&interaction"-field
++++++++++++++++++++


* **NORDER**-tag : The order of force constants to be calculated. Anharmonic terms up to :math:`(m+1)`\ th order will be considered with ``NORDER`` = :math:`m`.

 :Default: None
 :Type: Integer
 :Example: ``NORDER`` should be 1 for harmonic calculations, and 2 to include cubic terms.

````

* NBODY-tag : Entry for excluding multiple-body interactions from anharmonic force constants
 
 :Default: ``NBODY`` = [2, 3, 4, ..., ``NORDER`` + 1]
 :Type: Array of integers
 :Example: If one wants to exclude three-body interactions from cubic force constants, one should explicitly give ``NBODY = 2 2``.

````

"&cutoff"-field
+++++++++++++++

In this entry field, one needs to specify cutoff radii of interaction for each order in units of Bohr. 
In the current implementation, cutoff radii should be defined for every possible pair of atomic elements. 
For example, the cutoff entry for a harmonic calculation (``NORDER = 1``) of Si (``NKD = 1``) should be like
::

 &cutoff
  Si-Si 10.0
 /

This means that the cutoff radii of 10 :math:`a_{0}` will be used for harmonic Si-Si terms. 
Please note that the first column should be two character strings, which are contained in the ``KD``-tag, 
connected by a hyphen (’-’). 

When one wants to consider cubic terms (``NORDER = 2``), please specify the cutoff radius for cubic terms in the third column as the following::

 
 &cutoff
  Si-Si 10.0 5.6 # Pair r_{2} r_{3}
 /

Instead of giving specific cutoff radii, one can write "None" as follows::

 &cutoff
  Si-Si None 5.6
 /

which means that all possible harmonic terms between Si-Si atoms will be included. 

.. caution::

  Setting 'None' for anharmonic terms can greatly increase the number of parameters and thereby increase the computational cost.

When there are more than two atomic elements, please specify the cutoff radii between every possible pair of atomic elements. In the case of MgO (``NKD = 2``), the cutoff entry should be like
::
 
 &cutoff
  Mg-Mg 8.0
  O-O 8.0
  Mg-O 10.0
 /

which can equivalently be written by using the wild card (’*’) as
::

 &cutoff
  *-* 8.0
  Mg-O 10.0 # Overwrite the cutoff radius for Mg-O harmonic interactions
 /

.. important::

  Cutoff radii specified by an earlier entry will be overwritten by a new entry that comes later.

Once the cutoff radii are properly given, harmonic force constants
:math:`\Phi_{i,j}^{\mu,\nu}` satisfying :math:`r_{ij} \le r_{c}^{\mathrm{KD}[i]-\mathrm{KD}[j]}` will be searched.

In the case of cubic terms, force constants :math:`\Phi_{ijk}^{\mu\nu\lambda}` satisfying :math:`r_{ij} \le r_{c}^{\mathrm{KD}[i]-\mathrm{KD}[j]}`, :math:`r_{ik} \le r_{c}^{\mathrm{KD}[i]-\mathrm{KD}[k]}`, and
:math:`r_{jk} \le r_{c}^{\mathrm{KD}[j]-\mathrm{KD}[k]}` will be searched and determined by fitting.

````

"&cell"-field
+++++++++++++

Please give the cell parameters in this entry in units of Bohr as the following::

 &cell
  a
  a11 a12 a13
  a21 a22 a23
  a31 a32 a33
 /

The cell parameters are then given by :math:`\vec{a}_{1} = a \times (a_{11}, a_{12}, a_{13})`,
:math:`\vec{a}_{2} = a \times (a_{21}, a_{22}, a_{23})`, and :math:`\vec{a}_{3} = a \times (a_{31}, a_{32}, a_{33})`.

````

"&position"-field
+++++++++++++++++

In this field, one needs to specify the atomic element and fractional coordinate of atoms in the supercell. 
Each line should be
::

  ikd xf[1] xf[2] xf[3]

where `ikd` is an integer specifying the atomic element (`ikd` = 1, ..., ``NKD``) and `xf[i]` is the
fractional coordinate of an atom. There should be ``NAT`` such lines in the &position entry field.


````

"&fitting"-field
+++++++++++++++++

This field is necessary when ``MODE = fitting``.

* **DFILE**-tag : File name containing atomic displacements in Cartesian coordinate

 :Default: None
 :Type: String
 :Description: The format of ``DFILE`` can be found :ref:`here <label_format_DFILE>`

````

* **FFILE**-tag : File name containing atomic forces in Cartesian coordinate

 :Default: None
 :Type: String
 :Description: The format of ``FFILE`` can be found :ref:`here <label_format_DFILE>`

````

* **NDATA**-tag : Number of displacement-force data sets

 :Default: None
 :Type: Integer
 :Description: ``DFILE`` and ``FFILE`` should contain at least ``NDATA``:math:`\times` ``NAT`` lines.

````

* NSTART, NEND-tags : Specifies the range of data to be used for fitting

 :Default: ``NSTART = 1``, ``NEND = NDATA``
 :Type: Integer
 :Example: To use the data in the range of [20:30] out of 50 entries, the tags should be ``NSTART = 20`` and ``NEND = 30``.

````

* ICONST-tag = 0 | 1 | 2 | 3

 ===== =========================================================================
   0    No constraints
   1    Constraints for translational invariance will be imposed between IFCs.
   2   | In addition to ``ICONST = 1``, constraints for rotational invariance
       | will be imposed up to (``NORDER`` + 1)th order.
   3   | In addition to ``ICONST = 2``, constraints for rotational invariance
       | between (``NORDER`` + 1)th order and (``NORDER`` + 2)th order, which
       | are zero, will be considered. 
  11   | Same as ``ICONST = 1`` but the constraint is imposed algebraically
       | rather than numerically.
 ===== =========================================================================

 :Default: 1
 :Type: Integer
 :Description: See :ref:`this page<constraint_IFC>` for the numerical formulae.

````

* ROTAXIS-tag : Rotation axis used to estimate constraints for rotational invariance. This entry is necessary when ``ICONST = 2, 3``.

 :Default: None
 :Type: String
 :Example: When one wants to consider the rotational invariance around the :math:`x`\ -axis, one should give ``ROTAXIS = x``. If one needs additional constraints for the rotation around the :math:`y`\ -axis, ``ROTAXIS`` should be ``ROTAXIS = xy``. 

````

* FC2XML-tag : XML file to which the harmonic terms will be fixed upon fitting

 :Default: None
 :Type: String
 :Description: When ``FC2XML``-tag is given, harmonic force constants will be fixed to the values stored in the ``FC2XML`` file. This may be useful for optimizing cubic and higher-order terms without changing the harmonic terms. Please make sure that the number of harmonic terms in the new computational condition is the same as that in the ``FC2XML`` file.

````

* FC3XML-tag : XML file to which the cubic terms will be fixed upon fitting

 :Default: None
 :Type: String
 :Description: Same as the ``FC2XML``-tag, but ``FC3XML`` is to fix cubic force constants. 

````

.. _label_format_DFILE:

Format of DFILE and FFILE
~~~~~~~~~~~~~~~~~~~~~~~~~

The displacement-force data sets obtained by first-principles (or classical force-field) calculations
have to be saved to ``DFILE`` and ``FFILE`` to estimate IFCs with ``MODE = fitting``.
In ``DFILE``, please explicitly specify the atomic displacements :math:`u_{\alpha}(\ell\kappa)` **in units of Bohr** as follows:
 
.. math::
    :nowrap:
  
    \begin{eqnarray*}
     u_{x}(1) & u_{y}(1) & u_{z}(1) \\
     u_{x}(2) & u_{y}(2) & u_{z}(2) \\
     & \vdots & \\
     u_{x}(\mathrm{NAT}) & u_{y}(\mathrm{NAT}) & u_{z}(\mathrm{NAT})
    \end{eqnarray*}

When there are ``NAT`` atoms in the supercell and ``NDATA`` data sets, 
there should be  ``NAT`` :math:`\times` ``NDATA`` lines in the ``DFILE`` without blank lines.
In ``FFILE``, please specify the corresponding atomic forces **in units of Ryd/Bohr**.
