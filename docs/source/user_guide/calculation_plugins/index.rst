
.. _sec.quantumespresso:

Available calculation plugins
-----------------------------

Description
^^^^^^^^^^^
`Quantum Espresso`_ is a suite of open-source codes for electronic-structure calculations from first principles, based on density-functional theory, plane waves, and pseudopotentials, freely `available online`_.
Documentation of the code and its internal details can also be found in the `Quantum Espresso manual`_.

.. _Quantum Espresso: http://www.quantum-espresso.org/
.. _available online: http://qe-forge.org/gf/project/q-e/frs/?action=FrsReleaseBrowse&frs_package_id=18
.. _Quantum Espresso manual: http://www.quantum-espresso.org/users-manual/

Currently supported codes include:

* `pw.x`: Ground state properties, total energy, ionic relaxation, molecular dynamics, forces, etc...
* `cp.x`: Car-Parrinello molecular dynamics
* `ph.x`: Phonons from density functional perturbation theory
* `q2r.x`: Fourier transform the dynamical matrices in the real space
* `matdyn.x`: Fourier transform the dynamical matrices in the real space
* `neb.x`: Energy barriers and reaction pathways using the Nudged Elastic Band (NEB) method
* `dos.x`: Compute density of states (DOS)
* `projwfc.x`: Projects wavefunctions onto orthogonalized atomic wavefunctions
* `pw2wannier90.x`: Interface between Quantum ESPRESSO and the `Wannier90`_ code

.. _Wannier90: http://www.wannier.org/

Moreover, support for further codes can be implemented adapting the **namelist** plugin.

Plugins
^^^^^^^

.. toctree::
   :maxdepth: 4

   pw
   cp
   ph
   q2r
   matdyn
   neb
   dos
   projwfc
   pw2wannier90

