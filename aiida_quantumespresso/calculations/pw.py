# -*- coding: utf-8 -*-
from __future__ import absolute_import

import os
import six

from aiida import orm
from aiida.common.lang import classproperty
from aiida.plugins import factories
from aiida_quantumespresso.calculations import BasePwCpInputGenerator


class PwCalculation(BasePwCpInputGenerator):
    """`CalcJob` implementation for the pw.x code of Quantum ESPRESSO."""

    _automatic_namelists = {
        'scf': ['CONTROL', 'SYSTEM', 'ELECTRONS'],
        'nscf': ['CONTROL', 'SYSTEM', 'ELECTRONS'],
        'bands': ['CONTROL', 'SYSTEM', 'ELECTRONS'],
        'relax': ['CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS'],
        'md': ['CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS'],
        'vc-md': ['CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS', 'CELL'],
        'vc-relax': ['CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS', 'CELL'],
    }

    # Keywords that cannot be set by the user but will be set by the plugin
    _blocked_keywords = [
        ('CONTROL', 'pseudo_dir'),
        ('CONTROL', 'outdir'),
        ('CONTROL', 'prefix'),
        ('SYSTEM', 'ibrav'),
        ('SYSTEM', 'celldm'),
        ('SYSTEM', 'nat'),
        ('SYSTEM', 'ntyp'),
        ('SYSTEM', 'a'),
        ('SYSTEM', 'b'),
        ('SYSTEM', 'c'),
        ('SYSTEM', 'cosab'),
        ('SYSTEM', 'cosac'),
        ('SYSTEM', 'cosbc'),
    ]

    _use_kpoints = True

    # Default input and output files
    _DEFAULT_INPUT_FILE = 'aiida.in'
    _DEFAULT_OUTPUT_FILE = 'aiida.out'

    # Not using symlink in pw to allow multiple nscf to run on top of the same scf
    _default_symlink_usage = False

    @classproperty
    def xml_filepaths(cls):
        """Return a list of XML output filepaths relative to the remote working directory that should be retrieved."""
        filepaths = []

        for filename in cls.xml_filenames:
            filepath = os.path.join(cls._OUTPUT_SUBFOLDER, '{}.save'.format(cls._PREFIX), filename)
            filepaths.append(filepath)

        return filepaths

    @classmethod
    def define(cls, spec):
        super(PwCalculation, cls).define(spec)
        spec.input('metadata.options.input_filename', valid_type=six.string_types, default=cls._DEFAULT_INPUT_FILE)
        spec.input('metadata.options.output_filename', valid_type=six.string_types, default=cls._DEFAULT_OUTPUT_FILE)
        spec.input('metadata.options.parser_name', valid_type=six.string_types, default='quantumespresso.pw')
        spec.input('kpoints', valid_type=orm.KpointsData,
            help='kpoint mesh or kpoint path')
        spec.input('hubbard_file', valid_type=orm.SinglefileData, required=False,
            help='SinglefileData node containing the output Hubbard parameters from a HpCalculation')
        spec.output('output_parameters', valid_type=orm.Dict,
            help='The `output_parameters` output node of the successful calculation.')
        spec.output('output_structure', valid_type=orm.StructureData, required=False,
            help='The `output_structure` output node of the successful calculation if present.')
        spec.output('output_trajectory', valid_type=orm.TrajectoryData, required=False)
        spec.output('output_array', valid_type=orm.ArrayData, required=False,
            help='The `output_array` output node of the successful calculation if present.')
        spec.output('output_band', valid_type=orm.BandsData, required=False,
            help='The `output_band` output node of the successful calculation if present.')
        spec.output('output_kpoints', valid_type=orm.KpointsData, required=False)
        spec.output('output_atomic_occupations', valid_type=orm.Dict, required=False)
        spec.default_output_node = 'output_parameters'

        # Unrecoverable errors: resources like the retrieved folder or its expected contents are missing
        spec.exit_code(200, 'ERROR_NO_RETRIEVED_FOLDER',
            message='The retrieved folder data node could not be accessed.')
        spec.exit_code(201, 'ERROR_NO_RETRIEVED_TEMPORARY_FOLDER',
            message='The retrieved temporary folder could not be accessed.')
        spec.exit_code(210, 'ERROR_OUTPUT_STDOUT_MISSING',
            message='The retrieved folder did not contain the required stdout output file.')
        spec.exit_code(220, 'ERROR_OUTPUT_XML_MISSING',
            message='The retrieved folder did not contain the required required XML file.')
        spec.exit_code(221, 'ERROR_OUTPUT_XML_MULTIPLE',
            message='The retrieved folder contained multiple XML files.')

        # Unrecoverable errors: required retrieved files could not be read, parsed or are otherwise incomplete
        spec.exit_code(300, 'ERROR_OUTPUT_FILES',
            message='Both the stdout and XML output files could not be read or parsed.')
        spec.exit_code(310, 'ERROR_OUTPUT_STDOUT_READ',
            message='The stdout output file could not be read.')
        spec.exit_code(311, 'ERROR_OUTPUT_STDOUT_PARSE',
            message='The stdout output file could not be parsed.')
        spec.exit_code(312, 'ERROR_OUTPUT_STDOUT_INCOMPLETE',
            message='The stdout output file was incomplete.')
        spec.exit_code(320, 'ERROR_OUTPUT_XML_READ',
            message='The XML output file could not be read.')
        spec.exit_code(321, 'ERROR_OUTPUT_XML_PARSE',
            message='The XML output file could not be parsed.')
        spec.exit_code(322, 'ERROR_OUTPUT_XML_FORMAT',
            message='The XML output file has an unsupported format.')
        spec.exit_code(350, 'ERROR_UNEXPECTED_PARSER_EXCEPTION',
            message='The parser raised an unexpected exception.')

        # Significant errors but calculation can be used to restart
        spec.exit_code(400, 'ERROR_OUT_OF_WALLTIME',
            message='The calculation stopped prematurely because it ran out of walltime.')
        spec.exit_code(410, 'ERROR_ELECTRONIC_CONVERGENCE_NOT_REACHED',
            message='The electronic minimization cycle did not reach self-consistency.')
        spec.exit_code(500, 'ERROR_IONIC_CONVERGENCE_NOT_REACHED',
            message='The ionic minimization cycle did not converge for the given thresholds.')
        spec.exit_code(501, 'ERROR_IONIC_CONVERGENCE_REACHED_EXCEPT_IN_FINAL_SCF',
            message='Then ionic minimization cycle converged but the thresholds are exceeded in the final SCF.')
        spec.exit_code(502, 'ERROR_IONIC_CYCLE_EXCEEDED_NSTEP',
            message='The ionic minimization cycle did not converge after the maximum number of steps.')
        spec.exit_code(510, 'ERROR_IONIC_CYCLE_ELECTRONIC_CONVERGENCE_NOT_REACHED',
            message='The electronic minimization cycle failed during an ionic minimization cycle.')
        spec.exit_code(511, 'ERROR_IONIC_CONVERGENCE_REACHED_FINAL_SCF_FAILED',
            message='The ionic minimization cycle converged, but electronic convergence was not reached in the '
                    'final SCF.')
        spec.exit_code(520, 'ERROR_IONIC_CYCLE_BFGS_HISTORY_FAILURE',
            message='The ionic minimization cycle terminated prematurely because of two consecutive failures in the '
                    'BFGS algorithm.')
        spec.exit_code(521, 'ERROR_IONIC_CYCLE_BFGS_HISTORY_AND_FINAL_SCF_FAILURE',
            message='The ionic minimization cycle terminated prematurely because of two consecutive failures in the '
                    'BFGS algorithm and electronic convergence failed in the final SCF.')

    @classproperty
    def input_file_name_hubbard_file(cls):
        """
        The relative file name of the file containing the Hubbard parameters if they should
        be read from file instead of specified in the input file cards. Requires the
        aiida-quantumespresso-hp plugin to be installed
        """
        try:
            HpCalculation = factories.CalculationFactory('quantumespresso.hp')
        except Exception:
            raise RuntimeError('this is determined by the aiida-quantumespresso-hp plugin but it is not installed')

        return HpCalculation.input_file_name_hubbard_file

    @classmethod
    def input_helper(cls, *args, **kwargs):
        """
        Validate if the keywords are valid Quantum ESPRESSO pw.x keywords, and
        also helps in preparing the input parameter dictionary in a
        'standardized' form (e.g., converts ints to floats when required,
        or if the flag flat_mode is specified, puts the keywords in the right
        namelists).

        This function calls :py:func:`aiida_quantumespresso.calculations.helpers.pw_input_helper`,
        see its docstring for further information.
        """
        from . import helpers
        return helpers.pw_input_helper(*args, **kwargs)
