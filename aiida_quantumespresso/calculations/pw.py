# -*- coding: utf-8 -*-
"""
Plugin to create a Quantum Espresso pw.x file.
"""
# TODO: COPY OUTDIR FROM PREVIOUS CALCULATION! Should be an input node of type
# RemoteData (or maybe subclass it?).
# TODO: tests!
# TODO: DOC + implementation of SETTINGS
# TODO: preexec, postexec
# TODO: Check that no further parameters are passed in SETTINGS
# TODO: many cards missing: check and implement
#       e.g.: ['CONSTRAINTS', 'OCCUPATIONS']
# TODO: implement pre... and post... hooks to add arbitrary strings before
#       and after a namelist, and a 'final_string' (all optional); useful
#       for development when new cards are needed
# TODO: all a lot of logger.debug stuff
import os

from aiida.common.utils import classproperty
from aiida.orm import CalculationFactory
from aiida.orm.calculation.job import JobCalculation
from aiida.orm.data.array.kpoints import KpointsData
from aiida.orm.data.singlefile import SinglefileData
from aiida_quantumespresso.calculations import BasePwCpInputGenerator


class PwCalculation(BasePwCpInputGenerator, JobCalculation):
    """
    Main DFT code (PWscf, pw.x) of the Quantum ESPRESSO distribution.
    For more information, refer to http://www.quantum-espresso.org/
    """
    # false due to PWscf bug, could be set to true on versions >= 5.1.0
    _default_symlink_usage = False

    @classproperty
    def xml_filepaths(cls):
        """Returns a list of relative filepaths of XML files."""
        filepaths = []

        for filename in cls.xml_filenames:
            filepath = os.path.join(cls._OUTPUT_SUBFOLDER, '{}.save'.format(cls._PREFIX), filename)
            filepaths.append(filepath)
        
        print "xml_filepaths() is returning:", filepaths
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
            HpCalculation = CalculationFactory('quantumespresso.hp')
        except Exception as exception:
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
