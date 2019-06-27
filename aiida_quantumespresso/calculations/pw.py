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
        spec.output('output_parameters', valid_type=orm.Dict)
        spec.output('output_structure', valid_type=orm.StructureData, required=False)
        spec.output('output_trajectory', valid_type=orm.TrajectoryData, required=False)
        spec.output('output_array', valid_type=orm.ArrayData, required=False)
        spec.output('output_band', valid_type=orm.BandsData, required=False)
        spec.output('output_kpoints', valid_type=orm.KpointsData, required=False)
        spec.output('output_atomic_occupations', valid_type=orm.Dict, required=False)
        spec.default_output_node = 'output_parameters'

        spec.exit_code(
            100, 'ERROR_NO_RETRIEVED_FOLDER', message='The retrieved folder data node could not be accessed.')
        spec.exit_code(
            101, 'ERROR_NO_RETRIEVED_TEMPORARY_FOLDER', message='The retrieved temporary folder could not be accessed.')
        spec.exit_code(
            110, 'ERROR_READING_OUTPUT_FILE', message='The output file could not be read from the retrieved folder.')
        spec.exit_code(
            115, 'ERROR_MISSING_XML_FILE', message='The required XML file is not present in the retrieved folder.')
        spec.exit_code(
            116, 'ERROR_MULTIPLE_XML_FILES', message='The retrieved folder contains multiple XML files.')
        spec.exit_code(
            117, 'ERROR_READING_XML_FILE', message='The required XML file could not be read.')
        spec.exit_code(120, 'ERROR_INVALID_OUTPUT', message='The output file contains invalid output.')

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

        This function calls
        :py:func:`aiida_quantumespresso.calculations.helpers.pw_input_helper`,
        see its docstring for further information.
        """
        from . import helpers
        return helpers.pw_input_helper(*args, **kwargs)
