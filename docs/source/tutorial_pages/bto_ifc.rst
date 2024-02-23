.. _label_tutorial_sto_scph:

.. raw:: html

    <style> .red {color:red} </style>

.. role:: red

.. |Angstrom|   unicode:: U+00C5 

BaTiO\ :sub:`3` : Anharmonic interatomic force constants (IFCs)
---------------------------------------------------

This page explains how to calculate anharmonic interatomic force constants (IFCs) using ALAMODE, especially for strongly anharmonic materials.
The target materials is cubic BaTiO3\ :sub:`3`, which exhibits a strong lattice anharmonicity.

The example input files are provided in **example/BaTiO3/reference**.

Let's move to the example directory

.. code-block:: bash

  $ cd ${ALAMODE_ROOT}/example/BaTiO3


.. _tutorial_BTO_IFC_step1:

1. Generate the randomly displaced supercells
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We use the *ab initio* molecular dynamics (AIMD) calculations to generate the supercells with random displacements.

The example VASP inputs are provided in **example/BaTiO3/reference/1_vasp_md**.
In this tutorial, we assume that we have already performed the AIMD calculation and obtained :red:`vasprun.xml`.

Next, we generate the supercells with random atomic displacements from the AIMD trajectory.

First, we move to the **1_configurations** directory and copy :red:`vasprun.xml`

.. code-block:: bash

  $ cd 1_configurations
  $ cp ../1_vasp_md/vasprun.xml ./

:red:`POSCAR_ref_supercell` is the structure of the reference supercell for which we want to calculate the IFCs.

.. code-block:: bash
  $ python3 $TOOLSDIR/displace.py --VASP POSCAR_ref_supercell -md vasprun.xml -e 1001:5000:50 --random --mag 0.04 --prefix disp_aimd+random_

Here, the option ``-e 1001:5000:50`` means that we sample from the 1001-th snapshot to the 5000-th snapshot with the sampling-step of 50 time-steps.
Thus, 80 configurations are generated in this case.
The option ``--random --mag 0.04`` adds random displacements of around 0.04 |Angstrom| to each atom in the extracted snapshots to reduce correlations between successive snapshots.
``$TOOLSDIR`` is the path of the tools directory in the alamode.

.. note::

    We don't need high accuracy in the AIMD calculation because its purpose is just to generate the random structure.
    In fact, we use the following parameters, which are different from the subsequent calculation of the displacement-force data.
    
    * ENCUT = 400, EDIFF = 1.0E-6

    * 2x2x2 kmesh

    The temperature in the AIMD calculation is chosen so that the generated trajectory widely sample the 
    low-energy landscape of the potential energy surface. We choose 300 K in this case, which is comparable 
    or lower than the structural transition temperatures of the target material. 

.. note::

    The number of random configurations should be chosen so that the generated set of IFCs
    converges with respect to it.
    Ideally, one should check the convergence of the calculated physical quantities by changing
    the number of random configurations from which one extract the anharmonic IFCs.
    
    It depends on problems, but a rule of thumb tells us that 100~1000 configurations will do 
    for the calculation of cubic and quartic IFCs.
    One can reduce the number if one calculates only the cubic IFCs.

.. _tutorial_BTO_IFC_step2:

2. Generate the displacement-force data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We calculate the atomic forces for each random configuration generated in :ref:`step 1<tutorial_BTO_IFC_step1>`.

The other VASP input files (:red:`INCAR` and :red:`KPOINTS`) are provided in **example/BaTiO3/reference/2_vasp_dfset**.

After collecting the resultant :red:`vasprun.xml` of each calculations in **example/BaTiO3/reference/2_vasp_dfset**, 
generate the displacement-force data with the command

.. code-block::

  $ cd ${ALAMODE_ROOT}/example/BaTiO3
  $ cd 2_vasp_dfset
  $ cp ../1_configurations/POSCAR_ref_supercell ./
  $ python3 $TOOLSDIR/extract.py --VASP=POSCAR_ref_supercell vasprun*.xml > DFSET_AIMD_random

The generated :red:`DFSET_AIMD_random` stores the atomic displacements and the atomic forces in each configuration, 
from which we can calculate the anharmonic IFCs.

.. _tutorial_BTO_IFC_step3:

3. Cross validation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this tutorial, we assume that the harmonic force constants are already calculated. 
Please use the method explained in :ref:`here<label_tutorial_01>` for the calculation of harmonic IFCs.





