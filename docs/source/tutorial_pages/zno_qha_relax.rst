.. _label_tutorial_zno_qha_relax:

.. raw:: html

    <style> .red {color:red} </style>

.. role:: red

.. |Angstrom|   unicode:: U+00C5 

ZnO : A QHA-based structural optimization example
---------------------------------------------------

This page explains how to calculate crystal structures at finite temperatures(:math:`T`) based on the quasiharmonic approximation (QHA).

Let's move to the example directory

.. code-block:: bash

  $ cd ${ALAMODE_ROOT}/example/ZnO/qha_relax

.. _tutorial_ZnO_QHA_step1:

1. Prepare force constants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this tutorial, we assume that the harmonic and the anharmonic force constants are already calculated up to 4th order.
Here, please unzip all the xml files in **example/ZnO/qha_relax** and **example/ZnO/qha_relax/strain_IFC** before running the calculation.


.. note::
  We use the harmonic IFCs calculated in :math:`4\times 4\times 2` supercell 
  and the anharmonic IFCs calculated in :math:`3\times 3\times 2` supercell
  due to the large computational cost of the DFT calculations for the 
  anharmonic IFCs.
  This treatment is justified because the higher-order IFCs tend to be more localized 
  in real space.

.. _tutorial_ZnO_QHA_step2:

2. Prepare the additional input files.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We need to calculate the elastic constants, the strain-force coupling, and the order strain-harmonic-IFC coupling
to calculate the :math:`T`-dependence of the shape of the unit cell.
These input files must be placed in the directory named as the ``STRAIN_IFC_DIR``-tag, 
which is specified in ``&relax``-field in the input file of :red:`anphon`.

* The second-order elastic constants (SOEC) and the third-order elastic constants (TOEC) need to be calculated by yourself.
  The name of the input file must be :red:`elastic_constants.in`.
  The format of :red:`elastic_constants.in` is as follows.
  ::
    SOEC
    V_cell * C_xx,xx
    V_cell * C_xx,xy
    V_cell * C_xx,xz
    ...
    V_cell * C_zz,zz
    TOEC
    V_cell * C_xx,xx,xx
    V_cell * C_xx,xx,xy
    ...
    V_cell * C_zz,zz,zz

  :math:`V_{cell}` is the shape of the unit cell and 
  :math:`C_{\mu_1 \nu_1, \mu_2 \nu_2} = \frac{1}{V}\frac{\partial U}{\partial u_{\mu_1 \nu_1} \partial u_{\mu_2 \nu_2}}`,
  :math:`C_{\mu_1 \nu_1, \mu_2 \nu_2, \mu_3 \nu_3} = \frac{1}{V}\frac{\partial U}{\partial u_{\mu_1 \nu_1} \partial u_{\mu_2 \nu_2} \partial u_{\mu_3 \nu_3}}`
  are the elastic constants.
  The values in :red:`elastic_constants.in` are in Rydberg unit.

* The strain-force coupling can be calculated using the `strainIFCcoupling <https://github.com/r-masuki/strainIFCcoupling>`_ code.

  If the strain-force coupling is zero, i.e. the atomic force is zero when we apply finite strain with fixed fractional atomic coordinates, you can set ``RENORM_2TO1ST=0`` and abbreviate the corresponding input file.
  To use ``RENORM_2TO1ST=1``, you need to impose rotational invariance on the IFCs (See Appendix C of the `original paper <https://arxiv.org/abs/2302.04537>`_ for details of the proof), which is not recommended because it usually worsens the fitting error.

  Note that the input file of the strain-force coupling must be :red:`strain_force.in`.

  In this tutorial, the following block of the input file means that  
  if we apply the strain :math:`u_{xx}= 0.005`, the atomic forces that acts on the first atom is 
  :math:`(f_x, f_y, f_z) = (0.000000,  -0.034812,  -0.022224)` [Ry/Bohr], e.t.c.
  ::
    xx 0.005 1.0
    0.000000  -0.034812  -0.022224
    0.000000  0.034812  -0.022224
    0.000000  -0.039854  0.022224
    0.000000  0.039854  0.022224

  The meaning of the weight ``1.0`` is similar to that in the next paragraph.

  Note that if you use the `strainIFCcoupling <https://github.com/r-masuki/strainIFCcoupling>`_ code, you can obtain a set of input files that follows this format.

* The strain-harmonic-IFC coupling can be calculated using the `strainIFCcoupling <https://github.com/r-masuki/strainIFCcoupling>`_ code.
  The input files are :red:`strain_harmonic.in` and the related xml files.
  
  For example,
  ::
    xx 0.005 1.0 ZnO442_harmonic_xx_0005.xml

  in :red:`strain_harmonic.in` means that :red:`ZnO442_harmonic_xx_0005.xml` is the xml file of the harmonic IFCs with the strain :math:`u_{xx} = 0.005`. The weight is ``1.0`` in this case because we use one-sided difference to calculate 

  :math:`\frac{\partial \widetilde{\Phi}_{\mu\mu'}(0\alpha,R'\alpha')}{\partial u_{xx}} \simeq \frac{\widetilde{\Phi}(u_{xx} = 0.005) - \widetilde{\Phi}(u_{\mu \nu} = 0.0)}{0.005}`.
  
  If you use the central difference method 

  :math:`\frac{\partial \widetilde{\Phi}_{\mu\mu'}(0\alpha,R'\alpha')}{\partial u_{xx}} \simeq \frac{\widetilde{\Phi}(u_{xx} = 0.005) - \widetilde{\Phi}(u_{\mu \nu} = -0.005)}{0.005\times2}`

  :math:`= 0.5\times \frac{\widetilde{\Phi}(u_{xx} = 0.005) - \widetilde{\Phi}(u_{\mu \nu} = 0.0)}{0.005} + 0.5\times \frac{\widetilde{\Phi}(u_{xx} = -0.005) - \widetilde{\Phi}(u_{\mu \nu} = 0.0)}{-0.005}`,

  the corresponding :red:`strain_harmonic.in` would be like
  ::
    xx 0.005 0.5 ZnO442_harmonic_xx_0005.xml
    xx 0.005 0.5 ZnO442_harmonic_xx_minus_0005.xml

  with respective weights of ``0.5``. Note that :red:`ZnO442_harmonic_xx_minus_0005.xml` is not provided in this tutorial.

  For the off-diagonal strain,
  :: 
    yz 0.005 1.0 ZnO442_harmonic_yz_00025.xml
  
  means that :red:`ZnO442_harmonic_yz_00025.xml` is the set of harmonic IFCs with :math:`u_{yz} = u_{zy} = 0.005/2 = 0.0025`.

  Note that if you use the `strainIFCcoupling <https://github.com/r-masuki/strainIFCcoupling>`_ code, you can obtain a set of input files that follows this format.

.. _tutorial_ZnO_QHA_step3:

3. Prepare the input file.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The input file for the :red:`anphon` calclation is :red:`ZnO_qha_thermo.in`.

Run the calculation with 

.. code-block:: bash 

  $ ${ALAMODE_ROOT}/anphon/anphon ZnO_scph_thermo.in > ZnO_scph_thermo.log


.. _tutorial_ZnO_QHA_step4:

4. Analyze the calculation results.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We can plot the :math:`T`-dependence of the thermal strain, which is written in :red:`ZnO_qha.umn_tensor`, with 

.. code-block:: bash

  $ gnuplot plot.plt

to obtain the followin figure.

.. figure:: ../../img/ZnO_thermal_strain.png
  :scale: 30%
  :align: center

  The temperature-dependence of the thermal strain of ZnO. In this wurtzite case, :math:`u_{xx} = u_{yy} = a(T)/a(T=0)-1.0`, :math:`u_{zz} = c(T)/c(T=0)-1.0`, where :math:`a(T)` and :math:`c(T)` are the :math:`T`-dependent lengths of the :math:`a` and :math:`c`-axis respectively.

The ZSISA and the v-ZSISA results can be obtained by changing ``QHA_SCHEME``-tag in ``&qha``-field.

We can see that ZSISA accurately reproduces the :math:`T`-dependence of the shape of the unit cell.
v-ZSISA underestimates the anisotropy of the thermal expansion, while it gives a good estimation of the :math:`T`-dependence of the volume of the unit cell, which is consistent with the theorem proved in the `original paper <https://arxiv.org/abs/2302.04537>`_ Please see the paper for details of ZSISA and v-ZSISA.

We can also calculate the :math:`T`-induced change of the electric polarization by 

:math:`P_{\mu}(T) - P_{\mu}(T=0) =\frac{1}{V_{cell}} \sum_{\alpha \nu} Z^*_{\alpha \mu \nu} u^{(0)}_{\alpha \nu}+\sum_{\mu_1 \nu_1}d_{\mu, \mu_1 \nu_1} u_{\mu_1 \nu_1},`

where :math:`Z^*_{\alpha \mu \nu}` are the Born effective charges and :math:`d_{\mu, \mu_1 \nu_1}` are the piezoelectric tensors, which can be calcualted using DFPT in the reference structure. The :math:`T`-dependent atomic displacements :math:`u^{(0)}_{\alpha \nu}` and the strain tensor :math:`u_{\mu_1 \nu_1}` are written in :red:`ZnO_qha.atom_disp` and :red:`ZnO_qha.umn_tensor` respectively.
