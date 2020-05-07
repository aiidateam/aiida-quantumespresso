Workflows
=========

.. toctree::
   :maxdepth: 2

This section describes the different base workchains and utilities for workchain development.

Utilities
+++++++++

.. automodule:: aiida_quantumespresso.common.workchain.utils
   :members:

Base workchains
+++++++++++++++

.. currentmodule: aiida_quantumespresso.common.workchain.base.restart
.. autoclass:: aiida_quantumespresso.common.workchain.base.restart.BaseRestartWorkChain
   :members:

   .. autoattribute:: aiida_quantumespresso.common.workchain.base.restart.BaseRestartWorkChain._verbose
   .. autoattribute:: aiida_quantumespresso.common.workchain.base.restart.BaseRestartWorkChain._calculation_class
   .. autoattribute:: aiida_quantumespresso.common.workchain.base.restart.BaseRestartWorkChain._error_handler_entry_point


Full workchains
+++++++++++++++

.. currentmodule: aiida_quantumespresso.workflows
.. autoclass:: aiida_quantumespresso.workflows.ph.base.PhBaseWorkChain

.. autoclass:: aiida_quantumespresso.workflows.pw.base.PwBaseWorkChain

.. autoclass:: aiida_quantumespresso.workflows.pw.relax.PwRelaxWorkChain

.. autoclass:: aiida_quantumespresso.workflows.pw.bands.PwBandsWorkChain