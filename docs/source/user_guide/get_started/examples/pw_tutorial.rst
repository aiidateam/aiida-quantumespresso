.. _my-ref-to-pw-tutorial:

Running a PWscf calculation
===================================

.. toctree::
   :maxdepth: 2

This section will show how to launch a single PWscf (``pw.x`` executable) calculation.
It is assumed that you have already performed the installation, and that you already:

  - setup a computer (with ``verdi``);
  - installed Quantum ESPRESSO on your machine or a cluster;
  - setup the code and computer you want to use.
  - installed pseudo potentials (with ``aiida-pseudo install sssp``)

The classic ``pw.x`` input file
--------------------------------

This is the input file of Quantum ESPRESSO that we will try to execute.
It consists in the total energy calculation of a 5-atom cubic cell of BaTiO\ :sub:`3`\.
If the following input is not clear to you, please refer to the `Quantum ESPRESSO documentation <https://www.quantum-espresso.org/>`_.

::

    &CONTROL
      calculation = 'scf'
      outdir = './out/'
      prefix = 'aiida'
      pseudo_dir = './pseudo/'
      restart_mode = 'from_scratch'
      verbosity = 'high'
      wf_collect = .true.
    /
    &SYSTEM
      ecutrho =   2.4000000000d+02
      ecutwfc =   3.0000000000d+01
      ibrav = 0
      nat = 5
      ntyp = 3
    /
    &ELECTRONS
      conv_thr =   1.0000000000d-06
    /
    ATOMIC_SPECIES
    Ba     137.33    Ba.pbesol-spn-rrkjus_psl.0.2.3-tot-pslib030.UPF
    Ti     47.88     Ti.pbesol-spn-rrkjus_psl.0.2.3-tot-pslib030.UPF
    O      15.9994   O.pbesol-n-rrkjus_psl.0.1-tested-pslib030.UPF
    ATOMIC_POSITIONS angstrom
    Ba           0.0000000000       0.0000000000       0.0000000000
    Ti           2.0000000000       2.0000000000       2.0000000000
    O            2.0000000000       2.0000000000       0.0000000000
    O            2.0000000000       0.0000000000       2.0000000000
    O            0.0000000000       2.0000000000       2.0000000000
    K_POINTS automatic
    4 4 4 0 0 0
    CELL_PARAMETERS angstrom
          4.0000000000       0.0000000000       0.0000000000
          0.0000000000       4.0000000000       0.0000000000
          0.0000000000       0.0000000000       4.0000000000

Using the ``aiida-quantumespresso`` plugin, you can submit the same calculation via the following Python script:

.. literalinclude:: pw_short_example.py
    :language: python
    :start-after: start-marker

Not only will AiiDA track the provenance of the entire calculation, it will also take care of preparing the scheduler submission script, submitting the calculation on the cluster, and getting the results back when it's done.

In the following sections, we explain all aspects of this script step-by-step.

Running a pw calculation with AiiDA
-----------------------------------

Now we are going to prepare a script to submit a job to your local installation of AiiDA.
This example will be a rather long script: in fact it assumes that nothing in your database,
so that we will have to load everything, like the pseudopotential files and the structure.
In a more practical situation, you might load data from the database and
perform a small modification to re-use it.

Let's say that through the ``verdi`` command you have already installed
a cluster, say ``TheHive``, and that you also compiled Quantum ESPRESSO
on the cluster, and installed the code pw.x with ``verdi`` with label ``pw-6.3``
for instance, so that in the rest of this tutorial we will reference to the
code as ``pw-6.3@TheHive``.

Let's start writing the python script.
First of all, we need to load the configuration concerning your
particular installation, in particular, the details of your database installation::

  #!/usr/bin/env python
  from aiida import load_profile
  load_profile()

Code
----


Now we have to select the code. Note that in AiiDA the object 'code' in the
database is meant to represent a specific executable, i.e. a given
compiled version of a code. Every calculation in AiiDA is linked to a *code*, installed on
a specific *computer*.
This means that if you install Quantum ESPRESSO on two computers *A* and *B*,
you will need to have two different 'codes' in the database
(although the source of the code is the same, the binary file is different).

If you setup the code ``pw-6.3`` on machine ``TheHive`` correctly, then it is
sufficient to write::

  codename = 'pw-6.3@TheHive'
  from aiida.orm import Code
  code = Code.get_from_string(codename)

where in the last line we just load the database object representing the code.

.. note:: the ``.get_from_string()`` method is just a helper method for user
  convenience, but there are some weird cases that cannot be dealt in a
  simple way (duplicated labels, code names
  that are an integer number, code names containing the '@' symbol, ...): try
  to not do this! This is not an error, but does not allow to use the
  ``.get_from_string()`` method to get those calculations.
  You can use directly the ``.get()`` method, for instance::

    code = Code.get(label='pw-6.3', machinename='TheHive')

  or even more generally get the code from its (integer) ``PK``::

    code = load_node(PK)

Once you have a code, you can start to assemble the inputs to run
a PWscf calculation.

.. note:: To learn more about calculations and processes in AiiDA you can
  refer to the :ref:`aiida:topics:calculations` and :ref:`aiida:topics:processes` sections
  of the AiiDA manual
  Remember that what is shown here refers to the Quantum ESPRESSO plugin:
  different codes will in general required different inputs.


Preparing a calculation
------------------------

To make it easier to prepare the inputs of an AiiDA calculation, process classes
have a ``get_builder`` method, that simplify the access to the inputs of a calculation.

To use it, you can load the calculation class that you need through the
``CalculationFactory`` ::

  PwCalculation = CalculationFactory('quantumespresso.pw')
  builder = PwCalculation.get_builder()

In a similar way, you can use the ``get_builder`` utility from a code, with the advantage
that the *code* input gets automatically populated: ::

  builder = code.get_builder()


If you are using the builder to a verdi shell, you will
see all the inputs that are available for a calculation by pressing the ``TAB`` key
after typing ``builder.``. To understand what type of data is expected for a
particular input, you can append the ``?`` or the ``??`` at the end of an input. For example: ::

    >>> builder.structure?
    Type:        property
    String form: <property object at 0x7f58bdcc6728>
    Docstring:   {"name": "structure", "required": "True", "valid_type": "<class '`aiida.orm.nodes.data.structure.StructureData`'>", "help": "The input structure."}

In this case, the helper lets you know you what kind of data it expects
(:py:class:`StructureData <aiida.orm.nodes.data.structure.StructureData>`),
that it is a *required* input for the calculation, and a short string explaining what the *structure* input represents.

The plugin requires at least the presence of:

    - An input structure;
    - A k-points mesh, in the form of a :py:class:`KpointsData <aiida.orm.nodes.data.array.kpoints.KpointsData>` object;
    - Pseudopotential files for the atomic species present;
    - A parameters dictionary, that contains the details of the Quantum ESPRESSO calculation;


Other inputs are optional, for example:

    - *metadata* is a dictionary of inputs that modify slightly the behaviour of a processes;
    - *settings* is a :py:class:`Dict <aiida.orm.nodes.data.dict.Dict>` dictionary that provides access to more advanced, non-default feature of the code.




Structure
---------

We now proceed in setting up the structure.

.. note:: Here we discuss only the main features of structures in AiiDA, needed
    to run a Quantum ESPRESSO PW calculation.

    For more detailed information, give a look to the
    :ref:`aiida:topics:data_types:materials:structure`.

There are two ways to do that in AiiDA, a first one is to use the AiiDA Structure,
which we will explain in the following; the second choice is the
`Atomic Simulation Environment (ASE) <http://wiki.fysik.dtu.dk/ase/>`_
which provides excellent tools to manipulate structures (the ASE Atoms object
needs to be converted into an AiiDA Structure, see the note at the end of the section).

We first have to load the abstract object class that describes a structure.
We do it in the following way: we load the ``DataFactory``, which is a tool to load the classes
by their name, and then call ``StructureData`` the abstract class that we loaded.
Note that it is not yet a class instance!
If you are not familiar with the terminology of object programming,
you can look at `Wikipedia <http://en.wikipedia.org/wiki/Object_(computer_science)>`_ for
their short explanation: in common speech that one refers to *a* file as a class,
while *the* file is the object or the class instance. In other words,
the class is our definition of the object ``Structure``,
while its instance is what will be saved as an object in the database::

  from aiida.plugins import DataFactory
  StructureData = DataFactory('core.structure')

We define the cell with a 3x3 matrix (we choose the convention
where each ROW represents a lattice vector), which in this case is just a cube of size 4 Angstroms::

  alat = 4. # angstrom
  cell = [[alat, 0., 0.,],
          [0., alat, 0.,],
          [0., 0., alat,],
         ]

Now, we create the ``StructureData`` instance, assigning immediately the cell.
Then, we append to the empty crystal cell the atoms, specifying their element name and their positions::

  # BaTiO3 cubic structure
  s = StructureData(cell=cell)
  s.append_atom(position=(0.,0.,0.),symbols='Ba')
  s.append_atom(position=(alat/2.,alat/2.,alat/2.),symbols='Ti')
  s.append_atom(position=(alat/2.,alat/2.,0.),symbols='O')
  s.append_atom(position=(alat/2.,0.,alat/2.),symbols='O')
  s.append_atom(position=(0.,alat/2.,alat/2.),symbols='O')

To see more methods associated to the ``StructureData`` class,
look at the :ref:`aiida:topics:data_types:materials:structure` documentation on the AiiDA manual.

.. note:: When you create a node (in this case a ``StructureData`` node) as
  described above, you are just creating it in the computer memory, and not
  in the database. This is particularly useful to run tests without filling
  the AiiDA database with garbage.

  It is not necessary to store anything to the database at this point;
  if, however, you want to directly store the structure in the
  database for later use, you can just call the ``store()`` method of the Node::

    s.store()

For an extended tutorial about the creation of ``Structure`` objects,
check :ref:`this tutorial on the AiiDA-core documentation <aiida:topics:data_types:materials:structure>`.

.. note:: AiiDA also supports  ASE structures. Once you created your structure
  with ASE, in an object instance called say ``ase_s``, you can
  straightforwardly use it to create the AiiDA ``StructureData``, as::

    s = StructureData(ase=ase_s)

  and then save it ``s.store()``.


Parameters
----------

Now we need to provide also the parameters of a Quantum ESPRESSO calculation,
like the cutoff for the wavefunctions, some convergence threshold, and so on.
The Quantum ESPRESSO pw.x plugin requires to pass this information within a
``Dict`` object, that is a specific AiiDA data node that can store a
dictionary (even nested) of basic data types: integers, floats, strings, lists,
dates, ...
We first load the class through the ``DataFactory``,
just like we did for the ``Structure``.
Then we create the instance of the object ``parameter``.
To represent closely the structure of the Quantum ESPRESSO input file,
``Dict`` is a nested dictionary, at the first level the namelists
(capitalized), and then the variables with their values (in lower case).

Note also that numbers and booleans are written in Python, i.e. ``False`` and
not the Fortran string ``.false.``!
::

    Dict = DataFactory('core.dict')

    parameters = Dict({
        'CONTROL': {
            'calculation': 'scf',
            'restart_mode': 'from_scratch',
            'wf_collect': True,
        },
        'SYSTEM': {
            'ecutwfc': 30.,
            'ecutrho': 240.,
        },
        'ELECTRONS': {
            'conv_thr': 1.e-6,
        }
    })

.. note:: also in this case, we chose not to store the ``parameters`` node.
  If we wanted, we could even have done it in a single line::

    parameters = Dict({...}).store()

The experienced Quantum ESPRESSO user will have noticed also that a couple of variables
are missing: the *prefix*, the *pseudo directory* and the *scratch directory* are
reserved to the plugin, which will use default values, and there are specific
AiiDA methods to restart from a previous calculation.

Input parameters validation
///////////////////////////
The dictionary provided above is the standard format for storing the inputs
of Quantum ESPRESSO pw.x in the database. It is important to store the inputs
of different calculations in a consistent way because otherwise later querying
becomes impossible (e.g. if different units are used for the same flags,
if the same input is provided in different formats, and so on).

In the PW input plugin, we provide a function that will help you in
both validating the input, and creating the input in the expected format
without remembering in which namelists the keywords are located.

You can access this function as follows. First, you define the input dictionary::

    test_dict = {
        'CONTROL': {
            'calculation': 'scf',
        },
        'SYSTEM': {
            'ecutwfc': 30.,
        },
        'ELECTRONS': {
            'conv_thr': 1.e-6,
        }
    }

Then, you can verify if the input is correct by using the
:py:func:`~aiida_quantumespresso.calculations.helpers.pw_input_helper` function,
conveniently exposes also as a ``input_helper`` class method of the ``PwCalculation`` class::

  resdict = CalculationFactory('quantumespresso.pw').input_helper(test_dict, structure=s)

If the input is invalid, the function will raise a ``InputValidationError``
exception, and the error message will have a verbose explanation of the possible
errors, and in many cases it will suggest how to fix them. Otherwise, in ``resdict``
you will find the same dictionary you passed in input, potentially slightly modified
to fix some small mistakes (e.g., if you pass an integer value where a float is expected,
the type will be converted). You can then use the output for the input ``Dict`` node::

  parameters = Dict(resdict)

As an example, if you pass an incorrect input, e.g. the following where we have introduced
a few errors::

    test_dict = {
        'CONTROL': {
            'calculation': 'scf',
        },
        'SYSTEM': {
            'ecutwfc': 30.,
            'cosab': 10.,
            'nosym': 1,
        },
        'ELECTRONS': {
            'convthr': 1.e-6,
            'ecutrho': 100.
        }
    }

After running the ``input_helper`` method, you will get an exception with a message
similar to the following::

  QEInputValidationError: Errors! 4 issues found:
  * You should not provide explicitly keyword 'cosab'.
  * Problem parsing keyword convthr. Maybe you wanted to specify one of these: conv_thr, nconstr, forc_conv_thr?
  * Expected a boolean for keyword nosym, found <type 'int'> instead
  * Error, keyword 'ecutrho' specified in namelist 'ELECTRONS', but it should be instead in 'SYSTEM'

As you see, a quite large number of checks are done, and if a name is not provided, a list of
similar valid names is provided (e.g. for the wrong keyword "convthr" above).

There are a few additional options that are useful:

  - If you don't remember the namelists, you can pass a 'flat' dictionary, without
    namelists, and add the ``flat_mode=True`` option to ``input_helper``. Beside the usual
    validation, the function will reconstruct the correct dictionary to pass as input for
    the AiiDA Quantum ESPRESSO calculation. Example::

        test_dict_flat = {
            'calculation': 'scf',
            'ecutwfc': 30.,
            'conv_thr': 1.e-6,
        }
        resdict = CalculationFactory('quantumespresso.pw').input_helper(
            test_dict_flat, structure=s, flat_mode = True)

    and after running, ``resdict`` will contain::

        test_dict = {
            'CONTROL': {
                'calculation': 'scf',
            },
            'SYSTEM': {
                'ecutwfc': 30.,
            },
            'ELECTRONS': {
                'conv_thr': 1.e-6,
            }
        }

    where the namelists have been automatically generated.


  - You can pass a string with a specific version number for a feature that was added
    only in a given version. For instance::

     resdict = CalculationFactory('quantumespresso.pw').input_helper(
         test_dict, structure=s,version='5.3.0')

    If the specific version is not among those for which we have a list of valid parameters,
    the error message will tell you which versions are available.


.. note:: We will try to maintain the input_helper every time a new version of Quantum ESPRESSO
   is released, but consider the ``input_helper`` function as a utility, rather than the
   official way to provide the input -- the only officially supported way to provide
   an input to pw.x is through a direct dictionary, as described earlier in the section "Parameters".
   This applies in particular if you are using very old versions of Quantum ESPRESSO, or customized versions
   that accept different parameters.


Multi-dimensional variables
///////////////////////////
The input format of pw.x contains various keywords that do not simply take the format of a key value pair, but
rather there will some indices in the key itself. Take for example the ``Hubbard_U(i)`` keyword of the ``SYSTEM`` card.
The Hubbard U value needs to be applied to a specific species and therefore the index ``i`` is required to be able to
make this distinction. Note that the value of ``i`` needs to correspond to the index of the species to which the
Hubbard U value needs to be applied.

The ``PwCalculation`` plugin makes this easy as it will do the conversion from kind name to species index automatically.
This allows you to specify a Hubbard U value by using a dictionary notation, where the key is the kind name to which
it should be applied. For example, if you have a structure with the kind ``Co`` and what it to have a Hubbard U value,
one can add the following in the parameter data dictionary::

    parameters = {
        'SYSTEM': {
            'hubbard_u': {
                'Co': 4.5
            }
        }
    }

This part of the parameters dictionary will be transformed by the plugin into the following input file::

    &SYSTEM
        hubbard_u(1) = 4.5
    /
    ATOMIC_SPECIES
    Co     58.933195 Co_pbe_v1.2.uspp.F.UPF
    Li     6.941 li_pbe_v1.4.uspp.F.UPF
    O      15.9994 O_pbe_v1.2.uspp.F.UPF

Note that since ``Co`` is listed as the first atomic species, the index in the ``hubbard_u(1)`` keyword reflects this.
The usage of a dictionary where the keys correspond to a kind of the input structure, will work for any keyword where
the index should correspond to the index of the atomic species. Examples of keywords where this approach will work::

    angle1(i)
    angle2(i)
    hubbard_alpha(i)
    hubbard_beta(i)
    hubbard_j0(i)
    hubbard_u(i)
    london_c6(i)
    london_rvdw(i)
    starting_charge(i)
    starting_magnetization(i)

However, there are also keywords that require more than index, or where the single index actually does not correspond
to the index of an atomic species. The list of keywords that match this description::

    efield_cart(i)
    fixed_magnetization(i)
    hubbard_j(i,ityp)
    starting_ns_eigenvalue(m,ispin,I)

To allow one to define these keywords, one can use nested lists, where the first few elements constitute all the index
values and the final element corresponds to the actual value. For example the following::

    parameters = {
        'SYSTEM': {
            'starting_ns_eigenvalue': [
                [1, 1, 3, 3.5],
                [2, 1, 1, 2.8]
            ]
        }
    }

will result in the following input file::

    &SYSTEM
        starting_ns_eigenvalue(1,1,3) = 3.5
        starting_ns_eigenvalue(2,1,1) = 2.8
    /

Note that any of the values within the lists that correspond to a kind in the input structure, will be replaced with the
index of the corresponding atomic species. For example::

    hubbard_j: [
        [2, 'Ni', 3.5],
        [2, 'Fe', 7.4],
    ]

would be formatted as::

    hubbard_j(2, 1) = 3.5
    hubbard_j(2, 3) = 7.4

Assuming the input structure contained the kinds 'Ni' and 'Fe', which would have received the atomic species indices 1
and 3 in the ultimate input file, respectively.

.. note::

    Nota bene: The code will not verify that a keyword actually requires an atomic species index in a certain position,
    and will indiscriminately map the value to an atomic species index if that value corresponds to a kind name.


K-points mesh
-------------

The k-points have to be saved in another kind of data, namely ``KpointsData``::

  KpointsData = DataFactory('core.array.kpoints')
  kpoints = KpointsData()
  kpoints.set_kpoints_mesh([4,4,4])

In this case it generates a 4*4*4 mesh without offset. To add an offset one
can replace the last line by::

  kpoints.set_kpoints_mesh([4,4,4],offset=(0.5,0.5,0.5))

.. note:: Only offsets of 0 or 0.5 are possible (this is imposed by PWscf).

You can also specify kpoints manually, by inputing a list of points
in crystal coordinates (here they all have equal weights)::

    import numpy
    kpoints.set_kpoints([[i,i,0] for i in numpy.linspace(0,1,10)],
        weights = [1. for i in range(10)])


A Gamma point calculation can be submitted by providing the 'gamma_only' flag to
the options dictionary ::

   kpoints.set_kpoints_mesh([1,1,1])
   builder.settings = Dict({'gamma_only': True})

Pseudopotentials
----------------

There is still one missing piece of information, that is the
pseudopotential files, one for each element of the structure.

.. note::

   For a more extended documentation on how to import pseudopotentials in the database,
   and how to handle and instal pseudopotential families, you can find more
   information in :ref:`aiida:topics:data_types:materials:upf`.

The ``builder.pseudos`` input is a dictionary, where the keys are the names of
the elements, and the values are the ``UpfData`` objects stored in the database.

It is possible to specify manually which pseudopotential files to use
for each atomic species. However, for any practical use, it is convenient
to use the pseudopotential families.

If you got one installed, you can simply tell the calculation to use the
pseudopotential family with a given name, and AiiDA will take care of
linking the proper pseudopotentials to the calculation, one for each atomic
species present in the input structure. This can be done using::

  from aiida.orm.nodes.data.upf import get_pseudos_from_structure
  builder.pseudos = get_pseudos_from_structure(structure, <PSEUDOPOTENTIAL_FAMILY_NAME>)


.. note::

   The list of pseudopotential families installed in your database can be accessed
   by command line with ::

       verdi data upf listfamilies


Labels and comments
-------------------

Sometimes it is useful to attach some notes to the calculation,
that may help you later understand why you did such a calculation,
or note down what you understood out of it.
Comments are a special set of properties of the calculation, in the sense
that it is one of the few properties that can be changed, even after
the calculation has run.

These properties can be set in the ``metadata`` input of the calculation,
with a *label* (a short description) and *description* (longer) ::

    builder.metadata.label = 'My generic title'
    builder.metadata.description ' a PWscf calculation of BaTiO3'

.. note::
   The ``TAB`` expansion works also on nested properties: from ``builder.metadata.``
   you can explore the available options

Calculation resources
---------------------

General options that are independent on  the code or the plugin
are grouped under  ``builder.metadata.options``.
Here you can set up the resource that you want to allocate
to this calculation, that will be passed to the scheduler
that handles the queue on your computer: ::

    builder.metadata.options.resources = {'num_machines': 1}
    builder.metadata.options.max_wallclock_seconds = 1800

More options are available, and can be explored by expanding
``builder.metadata.options.`` + ``TAB``.




Launching the calculation
-------------------------


If we are satisfied with what you created, it is time to attach all the required inputs
to the calculation: ::

    builder.structure = structure
    builder.kpoints = kpoints
    builder.parameters = parameters


The data nodes do not need to be stored at this point, as AiiDA will store them upon submission.


To execute the calculation, there are two possible ways:

   - run: the calculation gets executed in the shell, locking it until it is finished;
   - submit: the calculation is handled to the AiiDA daemon, and it will be running in the background

The commands are explained more in detail in the :ref:`aiida:tutorial` documentation of AiiDA.

To run your calculation, you can execute: ::

    from aiida.engine import run
    results = run(builder)

where the ``results`` variable will contain  a dictionary containing all the nodes that were produced as output.
Alternatively, it is possible to use either ``run.get_node`` or ``run.get_pk`` methods to retrive more information
about the calculation node: ::


    from aiida.engine import run
    result, node = run.get_node(builder)
    result, pk = run.get_pk(builder)

where ``pk`` and ``node`` will contain, respectively, the ``node`` object of the calculation or its ``pk``.

To submit the calculation to the daemon, you can use instead ::

    from aiida.engine import submit
    calc = submit(builder)

Note that, in this case, ``calc`` is the calculation node, and *not* the result dictionary.

.. note::
    In order to inspect the inputs created by AiiDA without actually running the calculation,
    we can perform a *dry run* of the submission process: ::

        builder.metadata.dry_run = True
        builder.metadata.store_provenance = False

    This will create the input files, that are available for inspection.

.. note::
     You're not forced to assign the inputs through the builder: they can be provided as keywords argument
     when you launch the calculation, passing the calculation class as the first argument: ::

        run(PwCalculation, structure=s, pseudos=pseudos, kpoints = kpoints, ...)



The calculation results are directly accessible as a result of the ``run`` job, or they can be
easily accessed from ``calc.res.``, a shortcut to the ``calc.outputs.output_parameters`` dictionary: ::

    calc.res.energy
    calc.res.energy_units

to access the final energy and its units.



Summarizing, we created all the inputs needed by a PW calculation,
that are: parameters, kpoints, pseudopotential files and the structure.
We then created the calculation, where we specified that it is a PW calculation
and we specified the details of the remote cluster.


.. To know how to monitor and check the state of submitted calculations, check the
   :ref:<AiiDA-core documentation`aiida:calculation_state`>.


To continue the tutorial with the ``ph.x`` phonon code of Quantum ESPRESSO,
continue here: :ref:`my-ref-to-ph-tutorial`.




Script: source code
-------------------

In this section you'll find two scripts that do what explained in the tutorial.
The compact is a script with a minimal configuration required.
You can copy and paste it (or download it), modify the two strings ``codename``
and ``pseudo_family`` with the correct values, and execute it with::

  python pw_short_example.py

(It requires to have one family of pseudopotentials configured).

Download: :download:`this example script <pw_short_example.py>`

Advanced features
-----------------
For a list of advanced features that can be activated (change of the
command line parameters, blocking some coordinates, ...) you can refer
to :ref:`this section<pw-advanced-features>`
in the pw.x input plugin documentation.

Importing previously run Quantum ESPRESSO pw.x calculations: ``PwImmigrant``
----------------------------------------------------------------------------

Once you start using AiiDA to run simulations, we believe that you will find it
so convenient that you will use it for all your calculations.

At the beginning, however, you may have some calculations that you already have
run and are sitting in some folders, and that you want to import inside AiiDA.

This can be achieved with the ``PwImmigrant`` class described below,
for which you can find a tutorial :ref:`here <pwimmigrant-tutorial>`.
