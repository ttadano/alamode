.. _label_tutorial_sto_scph:

.. raw:: html

    <style> .red {color:red} </style>

.. role:: red

.. |Angstrom|   unicode:: U+00C5 

BaTiO\ :sub:`3` : An SCPH-based structural optimization example
---------------------------------------------------

This page explains how to calculate crystal structures at finite temperatures based on the SCPH theory.

The example input files are provided in **example/BaTiO3/reference**.

Let's move to the example directory

.. code-block:: bash

  $ cd ${ALAMODE_ROOT}/example/BaTiO3


.. _tutorial_BTO_scph_relax_step1:

1. Prepare force constants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this tutorial, we assume that the harmonic and the anharmonic force constants are already calculated up to 4th order.
Please see :ref:`the previous tutorial<label_tutorial_sto_scph>` for the calculation of anharmonic IFCs.

2. Prepare the input files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
