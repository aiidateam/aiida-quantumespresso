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
* **parent_calculation**, two parents calculations are required for EPW: a PW calculation the get the
  ground state wavefunction and a PH calculation to get the dynamical matrices and electron-phonon matrix elements

  Note: There are no direct links between calculations. The use_parent_calculation will set
  a link to the RemoteFolder attached to that calculation. Alternatively, the method **use_parent_folder**
  can be used to set this link directly.

* **kpoints**, class :py:class:`KpointsData <aiida.orm.nodes.data.array.kpoints.KpointsData>`
  Reciprocal space k-points on which to build the dynamical matrices and electron-phonon matrices.
  Has to be a homogeneous mesh.

* **qpoints**, class :py:class:`KpointsData <aiida.orm.nodes.data.array.kpoints.KpointsData>`
  Reciprocal space q-points on which to build the dynamical matrices and electron-phonon matrices.
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
      'INPUTEPW', 'iverbosity': file prefix
      'INPUTEPW', 'dvscf_dir': place where the dvscf from ph.x are.
      'INPUTEPW', 'amass': Atomic mass
      'INPUTEPW', 'nk1': k-mesh on b1
      'INPUTEPW', 'nk2': k-mesh on b2
      'INPUTEPW', 'nk3': k-mesh on b3
      'INPUTEPW', 'nq1': q-mesh on b1
      'INPUTEPW', 'nq2': q-mesh on b2
      'INPUTEPW', 'nq3': q-mesh on b3

* **settings**, class :py:class:`Dict <aiida.orm.nodes.data.dict.Dict>` (optional)
  An optional dictionary that activates non-default operations. Possible values are:

    *  **'PARENT_CALC_OUT_SUBFOLDER_NSCF'**: string. The subfolder of the parent non-self-consistent calculation
       scratch to be copied in the new scratch.
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
