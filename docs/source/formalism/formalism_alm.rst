Formalism for program alm
=========================

Interatomic force constants (IFCs)
----------------------------------

The starting point of the computational methodology is to approximate the potential energy of interacting atoms 
by a Taylor expansion with respect to atomic displacements by

.. math::
    :label: U_Taylor

    &U - U_{0} = \sum_{n=1}^{N} U_{n} = U_{1} + U_{2} + U_{3} + \cdots, \\
    &U_{n} = \frac{1}{n!} \sum_{\substack{\ell_{1}\kappa_{1}, \dots, \ell_{n}\kappa_{n} \\ \mu_{1},\dots, \mu_{n}}} \Phi_{\mu_{1}\dots\mu_{n}}(\ell_{1}\kappa_{1};\dots;\ell_{n}\kappa_{n}) \; u_{\mu_{1}}(\ell_{1}\kappa_{1})\cdots u_{\mu_{n}}(\ell_{n}\kappa_{n}).

Here, :math:`u_{\mu}(\ell\kappa)` is the atomic displacement of :math:`\kappa`\ th atom in the :math:`\ell`\ th unit cell along :math:`\mu`\ th direction, and :math:`\Phi_{\mu_{1}\dots\mu_{n}}(\ell_{1}\kappa_{1};\dots;\ell_{n}\kappa_{n})` is the :math:`n`\ th-order interatomic force constant (IFC).


Symmetry relationship between IFCs
----------------------------------

The are several relationships between IFCs which may be used to reduce the number of independence IFCs. 

* Permutation

 Firstly, IFCs should be invariant under the exchange of triplet :math:`(\ell,\kappa,\mu)`, e.g. 

 .. math::

  \Phi_{\mu_{1}\mu_{2}\mu_{3}}(\ell_{1}\kappa_{1};\ell_{2}\kappa_{2};\ell_{3}\kappa_{3})=\Phi_{\mu_{1}\mu_{3}\mu_{2}}(\ell_{1}\kappa_{1};\ell_{3}\kappa_{3};\ell_{2}\kappa_{2})=\dots. 

* Periodicity

 Secondly, since IFCs should depend on interatomic distances, they are invariant under a translation in units of lattice vector, namely 

 .. math::
 
  \Phi_{\mu_{1}\mu_{2}\dots\mu_{n}}(\ell_{1}\kappa_{1};\ell_{2}\kappa_{2};\dots;\ell_{n}\kappa_{n})=\Phi_{\mu_{1}\mu_{2}\dots\mu_{n}}(0\kappa_{1};\ell_{2}-\ell_{1}\kappa_{2};\dots;\ell_{n}-\ell_{1}\kappa_{n}). 

* Crystal symmetry

 A crystal symmetry operation maps an atom :math:`\vec{r}(\ell\kappa)` to another equivalent atom :math:`\vec{r}(LK)` by rotation and translation.
 Since the potential energy is invariant under any crystal symmetry operations, IFCs should transform under a symmetry operation as follows:

 .. math::
     :label: ifcsym1

     \sum_{\nu_{1},\dots,\nu_{n}}\Phi_{\nu_{1}\dots\nu_{n}}(L_{1}K_{1};\dots;L_{n}K_{n}) O_{\mu_{1}\nu_{1}}\cdots O_{\mu_{n}\nu_{n}} = \Phi_{\mu_{1}\dots\mu_{n}}(\ell_{1}\kappa_{1};\dots;\ell_{n}\kappa_{n}),

 where :math:`O` is the rotational matrix of the symmetry operation. 
 Let :math:`N_{s}` be the number of symmetry operations, there are :math:`N_{s}` relationships between IFCs which may be used to find independent IFCs.

 .. Note::

   In the current implementation of *alm*, independent IFCs are searched in Cartesian coordinate where the matrix element :math:`O_{\mu\nu}` is 0 or :math:`\pm1` in all symmetry operations except for those of **hexagonal** (trigonal) lattice. Also, except for hexagonal (trigonal) systems, the product :math:`O_{\mu_{1}\nu_{1}}\cdots O_{\mu_{n}\nu_{n}}` in the left hand side of equation :eq:`ifcsym1` becomes non-zero only for a specific pair of :math:`\{\nu\}` (and becomes 0 for all other :math:`\{\nu\}`\ s). Therefore, let :math:`\{\nu^{\prime}\}` be such a pair of :math:`\{\nu\}`, the equation :eq:`ifcsym1` can be reduced to

   .. math::
       :label: ifcsym2
     
       \Phi_{\nu_{1}^{\prime}\dots\nu_{n}^{\prime}}(L_{1}K_{1};\dots;L_{n}K_{n}) = s \Phi_{\mu_{1}\dots\mu_{n}}(\ell_{1}\kappa_{1};\dots;\ell_{n}\kappa_{n}),
   
   where :math:`s=\pm1`. The code employs equation :eq:`ifcsym2` instead of equation :eq:`ifcsym1` to reduce the number of IFCs. If IFCs of the left-hand side and the right-hand side of equation :eq:`ifcsym2` are equivalent and the coupling coefficient is :math:`s=-1`, the IFC is removed since it becomes zero. For **hexagonal** (trigonal) systems, there can be symmetry operations where multiple terms in the left-hand side of equation :eq:`ifcsym1` become non-zero. For such cases, equation :eq:`ifcsym1` is not used to reduce the number of IFCs. Alternatively, the corresponding symmetry relationships are imposed as constraints between IFCs in solving fitting problems.


.. _constraint_IFC:

Constraints between IFCs
------------------------

Since the potential energy is invariant under rigid translation and rotation, it may be necessary for IFCs to satisfy corresponding constraints.

The constraints for translational invariance are given by

.. math::
    :label: consttran

    \sum_{\ell_{1}\kappa_{1}}\Phi_{\mu_{1}\mu_{2}\dots\mu_{n}}(\ell_{1}\kappa_{1};\ell_{2}\kappa_{2};\dots;\ell_{n}\kappa_{n}) = 0,
  
which should be satisfied for arbitrary pairs of :math:`\ell_{2}\kappa_{2},\dots,\ell_{n}\kappa_{n}` and :math:`\mu_{1},\dots,\mu_{n}`. The code *alm* imposes equation :eq:`consttran` by default (``ICONST = 1``). 

The constraints for rotational invariance are

.. math::
    
    &\sum_{\ell^{\prime}\kappa^{\prime}} (\Phi_{\mu_{1}\dots\mu_{n}\nu}(\ell_{1}\kappa_{1};\dots;\ell_{n}\kappa_{n};\ell^{\prime}\kappa^{\prime}) r_{\mu}(\ell^{\prime}\kappa^{\prime}) 
    - \Phi_{\mu_{1}\dots\mu_{n}\mu}(\ell_{1}\kappa_{1};\dots;\ell_{n}\kappa_{n};\ell^{\prime}\kappa^{\prime}) r_{\nu}(\ell^{\prime}\kappa^{\prime})) \\
    & \hspace{10mm} + \sum_{\lambda = 1}^{n}\sum_{\mu_{\lambda}^{\prime}} \Phi_{\mu_{1}\dots\mu_{\lambda}^{\prime}\dots\mu_{n}}(\ell_{1}\kappa_{1};\dots;\ell_{\lambda}\kappa_{\lambda};\dots;\ell_{n}\kappa_{n}) (\delta_{\mu,\mu_{\lambda}}\delta_{\nu,\mu_{\lambda}^{\prime}} - \delta_{\nu,\mu_{\lambda}}\delta_{\mu,\mu_{\lambda}^{\prime}}) = 0,

which must be satisfied for arbitrary pairs of :math:`(\ell_{1}\kappa_{1},\dots,\ell_{n}\kappa_{n};\mu_{1},\dots,\mu_{n};\mu,\nu)`.
This is complicated since :math:`(n+1)`\ th-order IFCs (first line) are related to :math:`n`\ th-order IFCs (second line).

For example, the constraints for rotational invariance related to harmonic terms can be found as 

.. math::
    :label: constrot1

    &\sum_{\ell_{2}\kappa_{2}} (\Phi_{\mu_{1}\nu}(\ell_{1}\kappa_{1};\ell_{2}\kappa_{2})r_{\mu}(\ell_{2}\kappa_{2})-\Phi_{\mu_{1}\mu}(\ell_{1}\kappa_{1};\ell_{2}\kappa_{2})r_{\nu}(\ell_{2}\kappa_{2})) \notag \\
    & \hspace{10mm} + \Phi_{\nu}(\ell_{1}\kappa_{1})\delta_{\mu,\mu_{1}} - \Phi_{\mu}(\ell_{1}\kappa_{1})\delta_{\nu,\mu_{1}} = 0,

and

.. math::
    :label: constrot2 

    & \sum_{\ell_{3}\kappa_{3}} (\Phi_{\mu_{1}\mu_{2}\nu}(\ell_{1}\kappa_{1};\ell_{2}\kappa_{2};\ell_{3}\kappa_{3}) r_{\mu}(\ell_{3}\kappa_{3}) \notag
    - \Phi_{\mu_{1}\mu_{2}\mu}(\ell_{1}\kappa_{1};\ell_{2}\kappa_{2};\ell_{3}\kappa_{3}) r_{\nu}(\ell_{3}\kappa_{3})) \\
    & \hspace{10mm} 
    + \Phi_{\nu\mu_{2}}(\ell_{1}\kappa_{1};\ell_{2}\kappa_{2})\delta_{\mu,\mu_{1}} 
    - \Phi_{\mu\mu_{2}}(\ell_{1}\kappa_{1};\ell_{2}\kappa_{2})\delta_{\nu,\mu_{1}} \notag \\
    & \hspace{10mm} + \Phi_{\mu_{1}\nu}(\ell_{1}\kappa_{1};\ell_{2}\kappa_{2})\delta_{\mu,\mu_{2}}
    - \Phi_{\mu_{1}\mu}(\ell_{1}\kappa_{1};\ell_{2}\kappa_{2})\delta_{\nu,\mu_{2}} = 0.
  
When ``NORDER = 1``, equation :eq:`constrot1` will be considered if ``ICONST = 2``, whereas equation :eq:`constrot2` will be neglected. To further consider equation :eq:`constrot2`, please use ``ICONST = 3``, though it may enforce a number of harmonic IFCs to be zero since cubic terms don't exist in harmonic calculations (``NORDER = 1``).


.. _fitting_formalism:

Estimate IFCs by least-square fitting
-------------------------------------

The code **alm** extracts harmonic and anharmonic IFCs from a displacement-force data set by solving the following linear least-square problem:

.. math::
   :label: lsq

   \text{minimize} \ \ \chi^{2} = \sum_{t}^{m} \sum_{i} \|F_{i,t}^{\mathrm{DFT}} - F_{i,t}^{\mathrm{ALM}} \|^{2}.

Here, :math:`m` is the number of atomic configurations and the index :math:`i = (\ell,\kappa,\mu)` is the triplet of coordinates. 
The model force :math:`F_{i,t}^{\mathrm{ALM}}` is a linear function of IFCs :math:`\{\Phi\}` which can be obtained by differentiating :math:`U` (Eq. :eq:`U_Taylor`) by :math:`u_{i}`.
The parameters (IFCs) are determined so as to best mimic the atomic forces obtained by DFT calculations, :math:`F_{i,t}^{\mathrm{DFT}}`. 

To evaluate goodness of fit, **alm** reports the relative error :math:`\sigma` defined by

.. math::
   :label: fitting_error

   \sigma = \sqrt{\frac{\chi^{2}}{\sum_{t}^{m}\sum_{i} (F_{i,t}^{\mathrm{DFT}})^{2}}},

where the numerator is the residual of fit and the denominator is the square sum of DFT forces.

