.. _my-ref-to-ph-tutorial:

Phonon
======

.. toctree::
   :maxdepth: 2

In this chapter will get you through the launching of a phonon calculation with Quantum ESPRESSO, with ``ph.x``, a density functional perturbation theory code.
For this tutorial, it is required that you managed to launch the ``pw.x`` calculation, which is at the base of the phonon code; and of course it is assumed that you already know how to use the QE code.

The input of a phonon calculation can be very simple, but requires the output of a previous ``pw.x`` calculation.
Here we will try to compute the dynamical matrix on a mesh of points (actually consisting of a 1x1x1 mesh for brevity).
The input file that we want to create is this one::

    &INPUTPH
        fildyn = 'DYN_MAT/dynamical-matrix-'
        ldisp = .true.
        nq1 = 1
        nq2 = 1
        nq3 = 1
        outdir = './out/'
        prefix = 'aiida'
        tr2_ph =   1.0000000000d-08
        verbosity = 'high'
    /

Walkthrough
-----------

Setting up the input for a ``PhCalculation`` is typically simpler than doing so for the ``PwCalculation``.
The only novel thing you will have to learn is how to link the ``pw.x`` calculation as a a parent.
As before, we'll write the script step by step.

We first load a couple of useful modules, and load the AiiDA profile::

    from aiida import orm, load_profile
    load_profile()

We'll also assume you've managed to run a previous ``pw.x`` calculation with AiiDA, and know the PK of the corresponding calculation job.

You need to have compiled the code, either locally or on a remote machine, and have configured a new ``ph.x`` code for AiiDA in the same way as you set up the one for ``pw.x``.
Then we load the ``Code`` class-instance from the database::

    code_label = 'qe-7.0-ph@localhost'
    code = orm.load_code(code_label)

Here, the full label of the code we want to run is ``qe-7.0-ph@localhost``, indicating it is a code set up on the local ``localhost`` computer.
From the code instance, we can obtain a builder for the calculation using the ``get_builder()`` method::

    builder = code.get_builder()

Next, we'll populate this builder with the required inputs.

Parameters
^^^^^^^^^^

Just like the for the ``pw.x`` calculation, we have to provide the input parameters as an instance of a ``Dict`` node.
Again, ``Dict`` will simply represent a nested dictionary in the database, namelists at the first level, and then variables and values.
But this time of course, we need to use the variables of ``ph.x``.
For this example, the parameters will be exceedingly simple::

    parameters = orm.Dict(
        dict={
            'INPUTPH': {
                'tr2_ph' : 1.0e-8,
            }
        }
    )
    builder.parameters = parameters

You may be wondering why inputs like ``fildyn`` and ``prefix`` are not specified.
These are set automatically by the ``PhCalculation`` plugin, and hence do not need to be provided.
In fact, trying to run with these inputs in the parameters will fail, since they are blocked to only allow certain defaults and avoid issues with linking the ``ph.x`` calculation to e.g. future ``q2r.x`` calculations.
Similarly, the ``ldisp`` and ``nq1`` inputs are not specified through the ``parameters`` input, but rather by the ``qpoints`` one, as we'll show in the next section.

Q-points
^^^^^^^^

Next we'll specify the q-point mesh we want to use for the ``ph.x`` calculation.
As mentioned above, for this example we'll be running a single q-point (i.e. a 1x1x1 mesh)::

    q_points = orm.KpointsData()
    q_points.set_kpoints_mesh([1, 1, 1])

    builder.qpoints = q_points

The plugin will automatically convert this ``KpointsData`` input to the necessary tags in the input file, similar to what is done for the ``PwCalculation``.


Parent calculation
^^^^^^^^^^^^^^^^^^

The phonon calculation needs to start from an already completed ``pw.x`` do the perturbation theory calculation.
In order to do this, you will set the directory where you ran the ``pw.x`` calculation as a ``parent_folder`` of the ``ph.x`` calculation.

You first need to remember the PK of the parent calculation that you launched previously (let's say it's ``8394``).
Let's load this calculation job node and get the ``remote_folder`` where the calculation ran from its outputs::

    parent_calculation = orm.load_node(8394)
    parent_folder = parent_calculation.outputs.remote_folder

Next, we'll add this ``RemoteData`` node as an input to the ``builder``::

    builder.parent_folder = parent_folder

Now, the ``PhCalculation`` plugin will copy the required output of the ``pw.x`` calculation before starting the ``ph.x`` run.

Scheduler settings
^^^^^^^^^^^^^^^^^^

Finally, we'll set the walltime of our calculation, as well as the resources we want to use::

    builder.metadata.options.max_wallclock_seconds = 30 * 60 # 30 min
    builder.metadata.options.resources = {'num_machines': 1}

These can depend on whether the code has been set up locally or remotely, and the scheduler used, see :ref:`aiida:topics:schedulers`.

Execution
^^^^^^^^^

Now, everything is ready, and just like the ``PwCalculation``, you only need to instruct AiiDA to submit the calculation to the daemon::

    from aiida.engine import submit
    submit(builder)

Complete script
---------------

The tutorial above describes how to set up the basic inputs for a ``ph.x`` calculation step by step.
Below you can find a complete script for you to adapt and run directly.
Make sure not to forget to change the ``parent_pk`` and ``code_label``!

::

    from aiida import orm, load_profile
    load_profile()

    #####################
    # ADAPT TO YOUR NEEDS
    parent_pk = 8394
    code_label = 'qe-7.0-ph@localhost'
    #####################

    code = orm.load_code(code_label)

    parameters = orm.Dict(
        dict={
            'INPUTPH': {
                'tr2_ph' : 1.0e-8,
            }
        }
    )

    q_points = orm.KpointsData()
    q_points.set_kpoints_mesh([1, 1, 1])

    parent_calculation = orm.load_node(parent_pk)
    parent_folder = parent_calculation.outputs.remote_folder

    builder = code.get_builder()

    builder.parameters = parameters
    builder.qpoints = q_points
    builder.parent_folder = parent_folder
    builder.metadata.options.max_wallclock_seconds = 30 * 60 # 30 min
    builder.metadata.options.resources = {'num_machines': 1}

    from aiida.engine import submit
    submit(builder)
