# -*- coding: utf-8 -*-
"""
Tests for the ph input plugin.
"""

import os

import aiida
from aiida.common.exceptions import InputValidationError
from aiida.common.folders import SandboxFolder
from aiida.orm import CalculationFactory, DataFactory
from aiida.backends.testbase import AiidaTestCase
from aiida.orm import Code

PhCalc = CalculationFactory('quantumespresso.ph')
StructureData = DataFactory('structure')
ParameterData = DataFactory('parameter')
UpfData = DataFactory('upf')
KpointsData = DataFactory('array.kpoints')


class TestQEPHInputGeneration(AiidaTestCase):
    """
    Tests for QE ph.x
    """
    @classmethod
    def setUpClass(cls):
        super(TestQEPHInputGeneration, cls).setUpClass()
        cls.calc_params = {
            'computer': cls.computer,
            'resources': {
                'num_machines': 1,
                'num_mpiprocs_per_machine': 1
            }
        }

        cls.code = Code()
        cls.code.set_remote_computer_exec((cls.computer, '/x.x'))
        cls.code.store()

    def test_inputs(self):
        import logging

        parameters = ParameterData(dict={
            'INPUTPH': {
                'tr2_ph': 1.0e-8,
            }})

        qpoints = KpointsData()
        qpoints.set_kpoints_mesh([1, 1, 1])
        settings = None

        builder = PhCalc.get_builder()
        builder.code = self.code
        builder.qpoints = qpoints
        builder.settings = settings

        with SandboxFolder() as f:
            # I use the same SandboxFolder more than once because nothing
            # should be written for these failing tests

            # Missing required input nodes
            with self.assertRaises(InputValidationError):
                builder.submit_test(folder=f)
            builder.parameters = parameters

            builder.submit_test(folder=f)

