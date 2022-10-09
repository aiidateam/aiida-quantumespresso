===========
Get started
===========

.. _installation:requirements:

Requirements
============

To work with ``aiida-quantumespresso``, you should have:

* installed ``aiida-core``
* configured an AiiDA profile.

Please refer to the `documentation <https://aiida.readthedocs.io/projects/aiida-core/en/latest/intro/get_started.html>`_ of ``aiida-core`` for detailed instructions.


.. _installation:installation:

Installation
============

The Python package can be installed from the Python Package index `PyPI <https://pypi.org/>`_ or directly from the source:

.. tab-set::

    .. tab-item:: PyPI

        The recommended method of installation is to use the Python package manager |pip|_:

        .. code-block:: console

            $ pip install aiida-quantumespresso

        This will install the latest stable version that was released to PyPI.

    .. tab-item:: Source

        To install the package from source, first clone the repository and then install using |pip|_:

        .. code-block:: console

            $ git clone https://github.com/aiidateam/aiida-quantumespresso
            $ pip install -e aiida-quantumespresso

        The ``-e`` flag will install the package in editable mode, meaning that changes to the source code will be automatically picked up.


.. _installation:configuration:

Configuration
=============

To enable tab-completion for the command line interface, execute the following shell command (depending on the shell):

.. tab-set::

    .. tab-item:: bash

        .. code-block:: console

            $ eval "$(_AIIDA_QUANTUMESPRESSO_COMPLETE=bash_source aiida-quantumespresso)"

    .. tab-item:: zsh

        .. code-block:: console

            $ eval "$(_AIIDA_QUANTUMESPRESSO_COMPLETE=zsh_source aiida-quantumespresso)"

    .. tab-item:: fish

        .. code-block:: console

            $ eval (env _AIIDA_QUANTUMESPRESSO_COMPLETE=fish_source aiida-quantumespresso)


Place this command in your shell or virtual environment activation script to automatically enable tab completion when opening a new shell or activating an environment.
This file is shell specific, but likely one of the following:

* the startup file of your shell (``.bashrc``, ``.zsh``, ...), if aiida is installed system-wide
* the `activators <https://virtualenv.pypa.io/en/latest/user_guide.html#activators>`_ of your virtual environment
* a `startup file <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#saving-environment-variables>`_ for your conda environment

.. important::

    After having added the line to the start up script, make sure to restart the terminal or source the script for the changes to take effect.


.. _installation:setup:

Setup
=====

.. _installation:setup:computer:

Computer
--------

To run Quantum ESPRESSO calculations on a compute resource, the computer should first be set up in AiiDA.
This can be done from the command line interface (CLI) or the Python application programming interface (API).
In this example, we will set up the ``localhost``, the computer where AiiDA itself is running:

.. tab-set::

    .. tab-item:: CLI

        To set up a computer, use the ``verdi`` CLI of ``aiida-core``.

        .. code-block:: console

            $ verdi computer setup -n -L localhost -H localhost -T core.local -S core.direct

        After creating the localhost computer, configure it using:

        .. code-block:: console

            $ verdi computer configure core.local localhost -n --safe-interval 0

        Verify that the computer was properly setup by running:

        .. code-block:: console

            $ verdi computer test localhost


    .. tab-item:: API

        To setup a computer using the Python API, run the following code in a Python script or interactive shell:

        .. code-block:: python

            from aiida.orm import Computer

            computer = Computer(
                label='localhost',
                hostname='localhost',
                transport_type='core.local',
                scheduler_type='core.direct'
            ).store()
            computer.configure()

For more detailed information, please refer to the documentation `on setting up compute resources <https://aiida.readthedocs.io/projects/aiida-core/en/latest/howto/run_codes.html#how-to-set-up-a-computer>`_.


.. _installation:setup:code:

Code
----

To run a Quantum ESPRESSO code, it should first be setup in AiiDA.
This can be done from the command line interface (CLI) or the Python application programming interface (API).
In this example, we will setup the ``pw.x`` code that is installed on the computer where AiiDA is running:

.. tab-set::

    .. tab-item:: CLI

        To setup a particular Quantum ESPRESSO code, use the ``verdi`` CLI of ``aiida-core``.

        .. code-block:: console

            $ verdi code setup -n --on-computer -Y localhost -L pw -P quantumespresso.pw --remote-abs-path /path/to/pw.x

    .. tab-item:: API

        To setup particular Quantum ESPRESSO code using the Python API, run the following code in a Python script or interactive shell:

        .. code-block:: python

            from aiida.orm import Code

            computer = load_computer('localhost')
            code = Code(
                label='pw',
                remote_computer_exec=(computer, '/path/to/pw.x')
                input_plugin_name='quantumespresso.pw',
            ).store()

.. important::

    Make sure to replace ``/path/to/pw.x`` with the actual absolute path to the ``pw.x`` binary.

For more detailed information, please refer to the documentation `on setting up codes <https://aiida.readthedocs.io/projects/aiida-core/en/latest/howto/run_codes.html#how-to-setup-a-code>`_.


.. _installation:setup:pseudopotentials:

Pseudopotentials
----------------

Many Quantum ESPRESSO codes require pseudo potentials.
The simplest way of installing these is through the ``aiida-pseudo`` plugin package.
This should come as a dependency of ``aiida-quantumespresso`` and should already be installed.
If this is not the case, it can be installed using:

.. code-block:: console

    $ pip install aiida-pseudo

At a minimum, at least one pseudo potential family should be installed.
We recommend using the |SSSP|_:

.. code-block:: console

    $ aiida-pseudo install sssp

For more detailed information on installing other pseudo potential families, please refer to the documentation of |aiida-pseudo|_.


.. |pip| replace:: ``pip``
.. _pip: https://pip.pypa.io/en/stable/

.. |aiida-pseudo| replace:: ``aiida-pseudo``
.. _aiida-pseudo: https://aiida-pseudo.readthedocs.io/en/latest/index.html

.. |SSSP| replace:: Standard Solid-State Pseudopotentials (SSSP)
.. _SSSP: https://www.materialscloud.org/discover/sssp/table/efficiency
