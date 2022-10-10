
.. _howto:calculations:pw:

``pw.x``
--------

The ``pw.x`` code of Quantum ESPRESSO performs many different kinds of self-consistent calculations of electronic-structure properties within Density-Functional Theory (DFT),  using a plane-wave basis set and pseudopotentials.
Examples of these properties include ground-state energy and one-electron (Kohn-Sham) orbitals, atomic forces, stresses, and structural optimization, also with variable cell.

================== ===============================================================================
Plugin class       :class:`~aiida_quantumespresso.calculations.pw.PwCalculation`
Plugin entry point ``quantumespresso.pw``
================== ===============================================================================


How to launch a ``pw.x`` calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Below is a script with a basic example of how to run a ``pw.x`` calculation through the ``PwCalculation`` plugin that computes the electronic ground state of an fcc silicon crystal:

.. literalinclude:: ../../tutorials/include/scripts/run_pw_basic.py
    :language: python

Note that you may have to change the name of the code that is loaded using ``load_code`` and the pseudopotential family loaded with ``load_group``.


How to define input file parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``pw.x`` code supports many parameters that can be defined through the input file, as shown on the `official documentation <https://www.quantum-espresso.org/Doc/INPUT_PW.html>`_.
The parameters are divided into section or "cards".
Parameters that are part of cards that start with an ampersand (``&``) should be specified through the ``parameters`` input of the ``PwCalculation`` plugin.
The parameters are specified using a Python dictionary, where each card is its own sub-dictionary, for example:

.. code-block:: python

    parameters = {
        'CONTROL': {
            'calculation': 'scf'
        },
        'SYSTEM': {
            'smearing': 'gaussian'
        },
        'ELECTRONS': {
            'electron_maxstep': 10
        }
    }

The parameters dictionary should be wrapped in a :class:`~aiida.orm.nodes.data.dict.Dict` node and assigned to the ``parameters`` input of the process builder:

.. code-block:: python

    from aiida.orm import Dict, load_code
    builder = load_code('pw').get_builder()
    parameters = {
        ...
    }
    builder.parameters = Dict(parameters)

.. warning::

    There are a number of input parameters that *cannot* be set, as they will be automatically set by the plugin based on other inputs, such as the ``structure``.
    These include:

    * ``CONTROL.pseudo_dir``
    * ``CONTROL.outdir``
    * ``CONTROL.prefix``
    * ``SYSTEM.celldm``
    * ``SYSTEM.nat``
    * ``SYSTEM.ntyp``
    * ``SYSTEM.a``
    * ``SYSTEM.b``
    * ``SYSTEM.c``
    * ``SYSTEM.cosab``
    * ``SYSTEM.cosac``
    * ``SYSTEM.cosbc``

    Defining them anyway will result in an exception when launching the calculation.


Multidimensional parameters
............................

The input format of ``pw.x`` contains various keywords that do not simply take the format of a key value pair, but rather there will some indices in the key itself.
Take for example the |starting_magnetization|_ keyword of the ``SYSTEM`` card.
The starting magnetization value needs to be applied to a specific species and therefore the index ``i`` is required to be able to make this distinction.

The ``PwCalculation`` plugin makes this easy as it will do the conversion from kind name to species index automatically.
This allows you to specify a starting magnetization value by using a dictionary notation, where the key is the kind name to which it should be applied.
For example, if you have a structure with the kind ``Co`` and want it to have a given starting magnetization, one can add the following in the parameter data dictionary:

.. code-block:: python

    parameters = {
        'SYSTEM': {
            'starting_magnetization': {
                'Co': 4.5
            }
        }
    }

This part of the parameters dictionary will be transformed by the plugin into the following input file:

.. code-block::

    &SYSTEM
        starting_magnetization(1) = 4.5
    /
    ATOMIC_SPECIES
    Co     58.93 Co.UPF
    Li     6.941 Li.UPF
    O      15.99 O.UPF

Note that since ``Co`` is listed as the first atomic species, the index in the ``starting_magnetization(1)`` keyword reflects this.
The usage of a dictionary where the keys correspond to a kind of the input structure, will work for any keyword where the index should correspond to the index of the atomic species.
Examples of keywords where this approach will work are:

* ``angle1(i)``
* ``angle2(i)``
* ``hubbard_alpha(i)``
* ``hubbard_beta(i)``
* ``hubbard_j0(i)``
* ``hubbard_u(i)``
* ``london_c6(i)``
* ``london_rvdw(i)``
* ``starting_charge(i)``
* ``starting_magnetization(i)``

There are also keywords that require more than index, or where the single index actually does not correspond to the index of an atomic species, such as the |starting_ns_eigenvalue|_ parameters.
To allow one to define these keywords, one can use nested lists, where the first few elements constitute all the index values and the final element corresponds to the actual value.
For example the following:

.. code-block:: python

    parameters = {
        'SYSTEM': {
            'starting_ns_eigenvalue': [
                [1, 1, 3, 3.5],
                [2, 1, 1, 2.8]
            ]
        }
    }

will result in the following input file:

.. code-block::

    &SYSTEM
        starting_ns_eigenvalue(1,1,3) = 3.5
        starting_ns_eigenvalue(2,1,1) = 2.8
    /

Note that any of the values within the lists that correspond to a kind in the input structure, will be replaced with the index of the corresponding atomic species.
For example:

.. code-block:: python

    hubbard_j: [
        [2, 'Ni', 3.5],
        [2, 'Fe', 7.4],
    ]

would be formatted as:

.. code-block:: python

    hubbard_j(2, 1) = 3.5
    hubbard_j(2, 3) = 7.4

Assuming the input structure contained the kinds ``Ni`` and ``Fe``, which would have received the atomic species indices 1 and 3 in the ultimate input file, respectively.


How to define pseudopotentials
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each ``pw.x`` calculation requires a pseudopotential to be specified for each kind in the structure.
These pseudopotentials can be specified in the ``PwCalculation`` plugin through the ``pseudos`` input namespace.
This input takes a dictionary, where the keys are the kind names and the values are instances of the :class:`~aiida_pseudo.data.pseudo.upf.UpfData` data plugin of the |aiida-pseudo|_ plugin package.
For example, if the input structure is a ``GaAs`` crystal, the pseudopotentials could be specified as follows:

.. code-block:: python

    from aiida.orm import load_code
    from aiida_pseudo.data.pseudo.upf import UpfData

    # Assume we have a UPF for Ga and As on the local file system
    upf_ga = UpfData('/path/to/Ga.upf')
    upf_as = UpfData('/path/to/As.upf')

    builder = load_code('pw').get_builder()
    builder.pseudos = {
        'Ga': upf_ga,
        'As': upf_as,
    }

.. tip::

    We recommend using the pseudopotentials provided by the |SSSP|_.
    The |aiida-pseudo| package provides an easy and automated way to install them.
    Please refer to the :ref:`section on pseudopotential setup <installation:setup:pseudopotentials>` for details.

Getting pseudopotentials from a family
......................................

If pseudopotentials have been installed as a family using the |aiida-pseudo|_ package, they can be retrieved as follows:

.. code-block:: python

    from ase.build import bulk
    from aiida.orm import StructureData, load_code, load_group

    # Load the pseudopotential family whose pseudos to use
    family = load_group('SSSP/1.1/PBE/effiency')
    structure = StructureData(ase=bulk('GaAs', 'fcc', 5.4))

    builder = load_code('pw').get_builder()
    builder.pseudos = family.get_pseudos_from_structure(structure=structure)

The :meth:`~aiida_pseudo.groups.family.pseudo.PseudoPotentialFamily.get_pseudos` method will automatically return a dictionary with the pseudos necessary for the specified ``structure`` that can be immediately assigned to the ``pseudos`` input of the builder.


Getting recommended cutoffs from a family
.........................................

Certain pseudopotential families provided by |aiida-pseudo|_, such as the :class:`~aiida_pseudo.groups.family.sssp.SsspFamily`, provide recommended energy cutoffs for the wave function and charge density, for each pseudopotential they provide.
Similar to the pseudopotentials themselves, these can easily be retrieved for any given ``StructureData``:

.. code-block:: python

    from ase.build import bulk
    from aiida.orm import StructureData, load_code, load_group

    # Load the pseudopotential family whose pseudos to use
    family = load_group('SSSP/1.1/PBE/effiency')
    structure = StructureData(ase=bulk('GaAs', 'fcc', 5.4))

    builder = load_code('pw').get_builder()
    cutoff_wfc, cutoff_rho = family.get_recommended_cutoffs(structure=structure, unit='Ry')
    builder.parameters = {
        'SYSTEM': {
            'ecutwfc': cutoff_wfc,
            'ecutrho': cutoff_rho,
        }
    }

Be sure to specify the ``unit`` as ``Ry`` as that is the unit that ``pw.x`` will expect.


.. |starting_magnetization| replace:: ``starting_magnetization``
.. _starting_magnetization: https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm287

.. |starting_ns_eigenvalue| replace:: ``starting_ns_eigenvalue``
.. _starting_ns_eigenvalue: https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm466

.. |aiida-pseudo| replace:: ``aiida-pseudo``
.. _aiida-pseudo: https://aiida-pseudo.readthedocs.io/

.. |SSSP| replace:: Standard Solid-State Pseudopotentials (SSSP)
.. _SSSP: https://www.materialscloud.org/discover/sssp/table/efficiency
