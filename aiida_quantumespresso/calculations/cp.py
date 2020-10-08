# -*- coding: utf-8 -*-
"""Plugin to create a Quantum Espresso pw.x file.

TODO: COPY OUTDIR FROM PREVIOUS CALCULATION! Should be an input node of type
     RemoteData (or maybe subclass it?).
TODO: tests!
TODO: DOC + implementation of SETTINGS
TODO: preexec, postexec
TODO: Check that no further parameters are passed in SETTINGS
TODO: many cards missing: check and implement
      e.g.: ['CONSTRAINTS', 'OCCUPATIONS']
TODO: implement pre_... and post_... hooks to add arbitrary strings before
      and after a namelist, and a 'final_string' (all optional); useful
      for development when new cards are needed
TODO: all a lot of logger.debug stuff
"""
import os

from aiida import orm
from aiida.common.lang import classproperty

from aiida_quantumespresso.calculations import BasePwCpInputGenerator


class CpCalculation(BasePwCpInputGenerator):
    """`CalcJob` implementation for the cp.x code of Quantum ESPRESSO."""

    # Constants to use in the calculation
    _CP_READ_UNIT_NUMBER = 50
    _CP_WRITE_UNIT_NUMBER = 51
    _FILE_XML_PRINT_COUNTER_BASENAME = 'print_counter.xml'
    _FILE_XML_PRINT_COUNTER = os.path.join(
        BasePwCpInputGenerator._OUTPUT_SUBFOLDER,
        '{}_{}.save'.format(BasePwCpInputGenerator._PREFIX, _CP_WRITE_UNIT_NUMBER),
        _FILE_XML_PRINT_COUNTER_BASENAME,
    )

    # Input file "sections" that we are going to write by calculation type
    # The term namelist is part of FORTRAN's jargon
    _automatic_namelists = {
        'scf': ['CONTROL', 'SYSTEM', 'ELECTRONS'],
        'nscf': ['CONTROL', 'SYSTEM', 'ELECTRONS'],
        'relax': ['CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS'],
        'cp': ['CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS'],
        'vc-cp': ['CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS', 'CELL'],
        'vc-relax': ['CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS', 'CELL'],
        'vc-wf': ['CONTROL', 'SYSTEM', 'ELECTRONS', 'WANNIER'],
    }

    # Pieces of input that we won't allow users to set
    _blocked_keywords = [
        ('CONTROL', 'pseudo_dir'),  # set later
        ('CONTROL', 'outdir'),  # set later
        ('CONTROL', 'prefix'),  # set later
        ('SYSTEM', 'celldm'),
        ('SYSTEM', 'nat'),  # set later
        ('SYSTEM', 'ntyp'),  # set later
        ('SYSTEM', 'a'),
        ('SYSTEM', 'b'),
        ('SYSTEM', 'c'),
        ('SYSTEM', 'cosab'),
        ('SYSTEM', 'cosac'),
        ('SYSTEM', 'cosbc'),
        ('CONTROL', 'ndr', _CP_READ_UNIT_NUMBER),
        ('CONTROL', 'ndw', _CP_WRITE_UNIT_NUMBER),
    ]

    # In cp calculations we won't use kpoints data
    _use_kpoints = False

    # Use low verbosity for cp calculations
    _default_verbosity = 'low'

    _cp_ext_list = [
        'cel',
        'con',
        'eig',
        'evp',
        'for',
        'nos',
        'pol',
        'pos',
        'spr',
        'str',
        'the',
        'vel',
        'wfc',
    ]

    _internal_retrieve_list = [
        os.path.join(
            BasePwCpInputGenerator._OUTPUT_SUBFOLDER,
            '{}.{}'.format(BasePwCpInputGenerator._PREFIX, ext),
        ) for ext in _cp_ext_list
    ] + [_FILE_XML_PRINT_COUNTER]

    # in restarts, it will copy from the parent the following
    _restart_copy_from = os.path.join(
        BasePwCpInputGenerator._OUTPUT_SUBFOLDER,
        '{}_{}.save'.format(BasePwCpInputGenerator._PREFIX, _CP_WRITE_UNIT_NUMBER),
    )

    # in restarts, it will copy the previous folder in the following one
    _restart_copy_to = os.path.join(
        BasePwCpInputGenerator._OUTPUT_SUBFOLDER,
        '{}_{}.save'.format(BasePwCpInputGenerator._PREFIX, _CP_READ_UNIT_NUMBER),
    )

    @classproperty
    def xml_filepaths(cls):
        """Return a list of relative filepaths of XML files."""
        # pylint: disable=no-self-argument,not-an-iterable
        filepaths = []

        for filename in cls.xml_filenames:
            filepath = os.path.join(
                cls._OUTPUT_SUBFOLDER,
                '{}_{}.save'.format(cls._PREFIX, cls._CP_WRITE_UNIT_NUMBER),
                filename,
            )
            filepaths.append(filepath)

        return filepaths

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.input('metadata.options.parser_name', valid_type=str, default='quantumespresso.cp')
        spec.output('output_trajectory', valid_type=orm.TrajectoryData)
        spec.output('output_parameters', valid_type=orm.Dict)
        spec.default_output_node = 'output_parameters'

        spec.exit_code(300, 'ERROR_NO_RETRIEVED_FOLDER',
            message='The retrieved folder data node could not be accessed.')
        spec.exit_code(301, 'ERROR_NO_RETRIEVED_TEMPORARY_FOLDER',
            message='The retrieved temporary folder could not be accessed.')
        spec.exit_code(303, 'ERROR_MISSING_XML_FILE',
            message='The required XML file is not present in the retrieved folder.')
        spec.exit_code(304, 'ERROR_OUTPUT_XML_MULTIPLE',
            message='The retrieved folder contains multiple XML files.')
        spec.exit_code(310, 'ERROR_OUTPUT_STDOUT_READ',
            message='The stdout output file could not be read.')
        spec.exit_code(311, 'ERROR_OUTPUT_STDOUT_PARSE',
            message='The output file contains invalid output.')
        spec.exit_code(312, 'ERROR_OUTPUT_STDOUT_INCOMPLETE',
            message='The stdout output file was incomplete probably because the calculation got interrupted.')
        spec.exit_code(320, 'ERROR_OUTPUT_XML_READ',
            message='The required XML file could not be read.')
        spec.exit_code(330, 'ERROR_READING_POS_FILE',
            message='The required POS file could not be read.')
        spec.exit_code(340, 'ERROR_READING_TRAJECTORY_DATA',
            message='The required trajectory data could not be read.')
