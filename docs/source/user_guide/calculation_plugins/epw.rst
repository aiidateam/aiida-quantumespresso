EPW
++

Description
-----------
Plugin for the Quantum Espresso epw.x executable.

Supported codes
---------------
* tested from epw.x v5.2 onwards.

Inputs
------

* **parent_folder_nscf**, an PW calculation (nscf) is required for EPW to get the ground state wavefunction.

* **parent_folder_ph**, a PH calculation is required to get the dynamical matrices and electron-phonon matrix elements
  on the coarse grid.

* **kpoints**, class :py:class:`KpointsData <aiida.orm.nodes.data.array.kpoints.KpointsData>`
  Reciprocal space k-points on which to build the dynamical matrices and electron-phonon matrices.
  Has to be a homogeneous mesh.

* **qpoints**, class :py:class:`KpointsData <aiida.orm.nodes.data.array.kpoints.KpointsData>`
  Reciprocal space q-points on which to build the dynamical matrices and electron-phonon matrices.
  Has to be a homogeneous mesh.

* **kfpoints**, class :py:class:`KpointsData <aiida.orm.nodes.data.array.kpoints.KpointsData>`
  Reciprocal space fine k-points on which to interpolate the electron-phonon matrices.
  Has to be a homogeneous mesh.

* **qfpoints**, class :py:class:`KpointsData <aiida.orm.nodes.data.array.kpoints.KpointsData>`
  Reciprocal space fine q-points on which to interpolate the electron-phonon matrices.
  Has to be a homogeneous mesh.

* **parameters**, class :py:class:`Dict <aiida.orm.nodes.data.dict.Dict>`
  Input parameters of eph.x, as a nested dictionary, mapping the input of EPW.
  Example::

      {"INPUTEPW":{"elph": ".true."},
      }

  A full list of variables and their meaning is found in the `epw.x documentation`_.

  .. _epw.x documentation: http://epw.org.uk/Documentation/Inputs

  Following keywords are already taken care of by AiiDA::

      'INPUTEPW', 'outdir': scratch directory
      'INPUTEPW', 'prefix': file prefix
      'INPUTEPW', 'verbosity': file prefix
      'INPUTEPW', 'dvscf_dir': place where the dvscf from ph.x are.
      'INPUTEPW', 'amass': Atomic mass
      'INPUTEPW', 'nk1': coarse k-mesh on b1
      'INPUTEPW', 'nk2': coarse k-mesh on b2
      'INPUTEPW', 'nk3': coarse k-mesh on b3
      'INPUTEPW', 'nq1': coarse q-mesh on b1
      'INPUTEPW', 'nq2': coarse q-mesh on b2
      'INPUTEPW', 'nq3': coarse q-mesh on b3
      'INPUTEPW', 'nkf1': fine k-mesh on b1
      'INPUTEPW', 'nkf2': fine k-mesh on b2
      'INPUTEPW', 'nkf3': fine k-mesh on b3
      'INPUTEPW', 'nqf1': fine q-mesh on b1
      'INPUTEPW', 'nqf2': fine q-mesh on b2
      'INPUTEPW', 'nqf3': fine q-mesh on b3

* **settings**, class :py:class:`Dict <aiida.orm.nodes.data.dict.Dict>` (optional)
  An optional dictionary that activates non-default operations. Possible values are:

    *  **'NAMELISTS'**: list of strings. Specify all the list of Namelists to be
       printed in the input file.
    *  **'PARENT_FOLDER_SYMLINK'**: boolean # If True, create a symlnk to the scratch
       of the parent folder, otherwise the folder is copied (default: False)
    *  **'CMDLINE'**: list of strings. parameters to be put after the executable and before the input file.
       Example: ["-npool","4"] will produce `ph.x -npool 4 < aiida.in`

Outputs
-------
At present there is no parser for EPW produced files. Only the raw text output can be obtained.


Errors
------
No parsing of EPW is implemented at present.


Example run
-----------

#. Configure the ``pw.x``, ``ph.x`` and ``epw.x`` codes called here ``pw@localhost``, ``ph@localhost`` and ``epw@localhost``

#. Run scf calculation and retrieve the node PK (named here ``NODE_PK_SCF``)

::

  aiida-quantumespresso calculation launch pw -X pw@localhost -p   SSSP_1.1_efficiency

#. Run a phonon calculation and retrieve the node PK (named here ``NODE_PK_PH``)

::

  aiida-quantumespresso calculation launch ph -X ph@localhost  -C NODE_PK_SCF

#. Run an nscf calculation. This is not standard and needs to be done within the verdi shell:

::

  import os
  import numpy as np
  from aiida.engine import submit
  from aiida import orm

  PwCalculation = CalculationFactory('quantumespresso.pw')

  first_pw = load_node(2540)
  builder = first_pw.get_builder_restart()
  updated_parameters = builder.parameters.get_dict()
  updated_parameters['CONTROL']['calculation'] = 'nscf'
  updated_parameters['SYSTEM']['nbnd'] = 10

  KpointsData = DataFactory('core.array.kpoints')
  kpoints = KpointsData()

  klist = np.zeros((216, 3))
  tt = 0
  for ii in np.arange(0, 1, 1.0/6):
    for jj in np.arange(0, 1, 1.0/6):
      for kk in np.arange(0, 1, 1.0/6):
        klist[tt, :] = np.array([ii, jj, kk])
        tt += 1
  kpoints.set_kpoints(klist, cartesian = False, weights= np.ones(216)*1.0/216)
  kpoints.store()

  builder.kpoints = kpoints
  builder.parameters = Dict(updated_parameters)

  builder.parent_folder = first_pw.outputs.remote_folder

  submit(builder)

Record the PK number from that calculation (named here ``NODE_PK_NSCF``)

#. Run an EPW calculation

::

  aiida-quantumespresso calculation launch epw -X epw@localhost --pw-nscf-parent NODE_PK_NSCF  --ph-parent NODE_PK_PH

#. Retrive your data from the EPW calculation

::

  verdi process list -a
  verdi calcjob gotocomputer NODE_PK_EPW
