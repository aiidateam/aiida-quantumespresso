# -*- coding: utf-8 -*-
"""
Plugin to create a Quantum Espresso pw.x file.

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
from __future__ import absolute_import

import os

import six
from aiida import orm
from aiida.common.lang import classproperty
from aiida.engine import CalcJob

from aiida_quantumespresso.calculations import BasePwCpInputGenerator


class CpCalculation(BasePwCpInputGenerator, CalcJob):
    """`CalcJob` implementation for the cp.x code of Quantum ESPRESSO."""

    # Constants to use in the calculation
    _CP_READ_UNIT_NUMBER = 50
    _CP_WRITE_UNIT_NUMBER = 51
    _FILE_XML_PRINT_COUNTER_BASENAME = "print_counter.xml"
    _FILE_XML_PRINT_COUNTER = os.path.join(
        BasePwCpInputGenerator._OUTPUT_SUBFOLDER,
        "{}_{}.save".format(BasePwCpInputGenerator._PREFIX, _CP_WRITE_UNIT_NUMBER),
        _FILE_XML_PRINT_COUNTER_BASENAME,
    )

    # Input file "sections" that we are going to write by calculation type
    # The term namelist is part of FORTRAN's jargon
    _automatic_namelists = {
        'scf': ["CONTROL", "SYSTEM", "ELECTRONS"],
        'nscf': ["CONTROL", "SYSTEM", "ELECTRONS"],
        'relax': ["CONTROL", "SYSTEM", "ELECTRONS", "IONS"],
        'cp': ["CONTROL", "SYSTEM", "ELECTRONS", "IONS"],
        'vc-cp': ["CONTROL", "SYSTEM", "ELECTRONS", "IONS", "CELL"],
        'vc-relax': ["CONTROL", "SYSTEM", "ELECTRONS", "IONS", "CELL"],
        'vc-wf': ["CONTROL", "SYSTEM", "ELECTRONS", "WANNIER"],
    }

    # Pieces of input that we won't allow users to set
    _blocked_keywords = [
        ("CONTROL", "pseudo_dir"),  # set later
        ("CONTROL", "outdir"),  # set later
        ("CONTROL", "prefix"),  # set later
        ("SYSTEM", "ibrav"),  # set later
        ("SYSTEM", "celldm"),
        ("SYSTEM", "nat"),  # set later
        ("SYSTEM", "ntyp"),  # set later
        ("SYSTEM", "a"),
        ("SYSTEM", "b"),
        ("SYSTEM", "c"),
        ("SYSTEM", "cosab"),
        ("SYSTEM", "cosac"),
        ("SYSTEM", "cosbc"),
        ("CONTROL", "ndr", _CP_READ_UNIT_NUMBER),
        ("CONTROL", "ndw", _CP_WRITE_UNIT_NUMBER),
    ]

    # In cp calculations we won't use kpoints data
    _use_kpoints = False

    # Use low verbosity for cp calculations
    _default_verbosity = "low"

    _cp_ext_list = [
        "cel",
        "con",
        "eig",
        "evp",
        "for",
        "nos",
        "pol",
        "pos",
        "spr",
        "str",
        "the",
        "vel",
        "wfc",
    ]

    _internal_retrieve_list = [
        os.path.join(
            BasePwCpInputGenerator._OUTPUT_SUBFOLDER,
            "{}.{}".format(BasePwCpInputGenerator._PREFIX, ext),
        )
        for ext in _cp_ext_list
    ] + [_FILE_XML_PRINT_COUNTER]

    # in restarts, it will copy from the parent the following
    _restart_copy_from = os.path.join(
        BasePwCpInputGenerator._OUTPUT_SUBFOLDER,
        "{}_{}.save".format(BasePwCpInputGenerator._PREFIX, _CP_WRITE_UNIT_NUMBER),
    )

    # in restarts, it will copy the previous folder in the following one
    _restart_copy_to = os.path.join(
        BasePwCpInputGenerator._OUTPUT_SUBFOLDER,
        "{}_{}.save".format(BasePwCpInputGenerator._PREFIX, _CP_READ_UNIT_NUMBER),
    )

    # Default options
    _DEFAULT_PARSER = "quantumespresso.cp"
    _DEFAULT_INPUT_FILE = "aiida.in"
    _DEFAULT_OUTPUT_FILE = "aiida.out"

    @classproperty
    def xml_filepaths(cls):
        """Returns a list of relative filepaths of XML files."""
        filepaths = []

        for filename in cls.xml_filenames:
            filepath = os.path.join(
                cls._OUTPUT_SUBFOLDER,
                "{}_{}.save".format(cls._PREFIX, cls._CP_WRITE_UNIT_NUMBER),
                filename,
            )
            filepaths.append(filepath)

        return filepaths

    @classmethod
    def define(cls, spec):
        super(CpCalculation, cls).define(spec)
        spec.input(
            "metadata.options.input_filename",
            valid_type=six.string_types,
            default=cls._DEFAULT_INPUT_FILE,
            non_db=True,
        )
        spec.input(
            "metadata.options.output_filename",
            valid_type=six.string_types,
            default=cls._DEFAULT_OUTPUT_FILE,
            non_db=True,
        )
        spec.input(
            "metadata.options.parser_name",
            valid_type=six.string_types,
            default=cls._DEFAULT_PARSER,
            non_db=True,
        )
        spec.output("output_trajectory", valid_type=orm.TrajectoryData)
        spec.output("output_parameters", valid_type=orm.Dict)
        spec.default_output_node = "output_parameters"

        spec.exit_code(
            100,
            "ERROR_NO_RETRIEVED_FOLDER",
            message="The retrieved folder data node could not be accessed.",
        )
        spec.exit_code(
            101,
            "ERROR_NO_RETRIEVED_TEMPORARY_FOLDER",
            message="The retrieved temporary folder could not be accessed.",
        )
        spec.exit_code(
            110,
            "ERROR_READING_OUTPUT_FILE",
            message="The output file could not be read from the retrieved folder.",
        )
        spec.exit_code(
            115,
            "ERROR_MISSING_XML_FILE",
            message="The required XML file is not present in the retrieved folder.",
        )
        spec.exit_code(
            116,
            "ERROR_MULTIPLE_XML_FILES",
            message="The retrieved folder contains multiple XML files.",
        )
        spec.exit_code(
            117,
            "ERROR_READING_XML_FILE",
            message="The required XML file could not be read.",
        )
        spec.exit_code(
            118,
            "ERROR_READING_POS_FILE",
            message="The required POS file could not be read.",
        )
        spec.exit_code(
            119,
            "ERROR_READING_TRAJECTORY_DATA",
            message="The required trajectory data could not be read.",
        )
        spec.exit_code(
            120,
            "ERROR_INVALID_OUTPUT",
            message="The output file contains invalid output.",
        )
