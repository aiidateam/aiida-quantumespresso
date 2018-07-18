# -*- coding: utf-8 -*-
from __future__ import absolute_import
import os

import six

from aiida import orm
from aiida.common.lang import classproperty
from aiida.plugins import factories
from aiida_quantumespresso.calculations import BasePwCpInputGenerator


class PwCalculation(BasePwCpInputGenerator):

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

    _DATAFILE_XML = os.path.join(
        BasePwCpInputGenerator._OUTPUT_SUBFOLDER,
        '{}.save'.format(BasePwCpInputGenerator._PREFIX),
        BasePwCpInputGenerator._DATAFILE_XML_BASENAME)

    _use_kpoints = True

    # Default input and output files
    _DEFAULT_INPUT_FILE = 'aiida.in'
    _DEFAULT_OUTPUT_FILE = 'aiida.out'

    #not using symlink in pw to allow multiple nscf to run on top of the same scf
    _default_symlink_usage = False

    @classmethod
    def define(cls, spec):
        super(PwCalculation, cls).define(spec)
        spec.input('metadata.options.input_filename', valid_type=six.string_types, default=cls._DEFAULT_INPUT_FILE, non_db=True)
        spec.input('metadata.options.output_filename', valid_type=six.string_types, default=cls._DEFAULT_OUTPUT_FILE, non_db=True)
        spec.input('metadata.options.parser_name', valid_type=six.string_types, default='quantumespresso.pw', non_db=True)
        spec.input('kpoints', valid_type=orm.KpointsData,
            help='kpoint mesh or kpoint path')
        spec.input('hubbard_file', valid_type=orm.SinglefileData, required=False,
            help='SinglefileData node containing the output Hubbard parameters from a HpCalculation')
        spec.exit_code(
            100, 'ERROR_NO_RETRIEVED_FOLDER', message='The retrieved folder data node could not be accessed.')
        spec.exit_code(
            110, 'ERROR_READING_OUTPUT_FILE', message='The output file could not be read from the retrieved folder.')
        spec.exit_code(120, 'ERROR_INVALID_OUTPUT', message='The output file contains invalid output.')

    @classproperty
    def xml_filepaths(cls):
        """Returns a list of relative filepaths of XML files."""
        filepaths = []

        for filename in cls.xml_filenames:
            filepath = os.path.join(cls._OUTPUT_SUBFOLDER, '{}.save'.format(cls._PREFIX), filename)
            filepaths.append(filepath)

        return filepaths

    def _init_internal_params(self):
        super(PwCalculation, self)._init_internal_params()

        self._xml_files = []

        # Default PW output parser provided by AiiDA
        self._default_parser = 'quantumespresso.pw'

        self._automatic_namelists = {
            'scf': ['CONTROL', 'SYSTEM', 'ELECTRONS'],
            'nscf': ['CONTROL', 'SYSTEM', 'ELECTRONS'],
            'bands': ['CONTROL', 'SYSTEM', 'ELECTRONS'],
            'relax': ['CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS'],
            'md': ['CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS'],
            'vc-md': ['CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS', 'CELL'],
            'vc-relax': ['CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS', 'CELL'],
        }

        # Keywords that cannot be set by the user but will be set by the plugin
        self._blocked_keywords = [
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

        self._use_kpoints = True

        # Default input and output files
        self._DEFAULT_INPUT_FILE = 'aiida.in'
        self._DEFAULT_OUTPUT_FILE = 'aiida.out'

    @classproperty
    def _use_methods(cls):
        """
        Extend the parent _use_methods with further keys.
        """
        use_methods = JobCalculation._use_methods
        use_methods.update(BasePwCpInputGenerator._baseclass_use_methods)
        use_methods.update({
            'kpoints': {
                'valid_types': (KpointsData, ),
                'additional_parameter': None,
                'linkname': 'kpoints',
                'docstring': ('KpointsData node that contains the kpoints mesh or path'),
            },
            'hubbard_file': {
                'valid_types': (SinglefileData, ),
                'additional_parameter': None,
                'linkname': 'hubbard_file',
                'docstring': 'SinglefileData node containing the output Hubbard parameters from a HpCalculation',
            }
        })

        return use_methods

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
