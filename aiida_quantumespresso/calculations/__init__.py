# -*- coding: utf-8 -*-
from __future__ import absolute_import
import io
import abc
import collections
import os

from aiida import orm
from aiida.common import datastructures, exceptions
from aiida.engine import CalcJob
from aiida.common.exceptions import InputValidationError
from aiida.common.datastructures import CalcInfo
from aiida.orm.data.upf import get_pseudos_from_structure
from aiida.common.utils import classproperty
from aiida.orm.data.structure import StructureData
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.array.kpoints import KpointsData
from aiida.orm.data.upf import UpfData
from aiida.orm.data.singlefile import SinglefileData
from aiida.orm.data.remote import RemoteData
from aiida.common.datastructures import CodeInfo
from aiida_quantumespresso.utils.convert import convert_input_to_namelist_entry

import six
from six.moves import zip

class BasePwCpInputGenerator(CalcJob):

    _PSEUDO_SUBFOLDER = './pseudo/'
    _OUTPUT_SUBFOLDER = './out/'
    _PREFIX = 'aiida'
    _INPUT_FILE_NAME = 'aiida.in'
    _OUTPUT_FILE_NAME = 'aiida.out'
    _DATAFILE_XML_PRE_6_2 = 'data-file.xml'
    _DATAFILE_XML_POST_6_2 = 'data-file-schema.xml'
    _ENVIRON_INPUT_FILE_NAME = 'environ.in'

    # Additional files that should always be retrieved for the specific plugin
    _internal_retrieve_list = []

    # Default PW output parser provided by AiiDA to be defined in the subclass
    _automatic_namelists = {}

    # Blocked keywords that are to be specified in the subclass
    _blocked_keywords = {}

    # In restarts, will not copy but use symlinks
    _default_symlink_usage = True

    # In restarts, it will copy from the parent the following
    _restart_copy_from = os.path.join(_OUTPUT_SUBFOLDER, '*')

    # In restarts, it will copy the previous folder in the following one
    _restart_copy_to = _OUTPUT_SUBFOLDER

    # Default verbosity; change in subclasses
    _default_verbosity = 'high'

    @classproperty
    def xml_filenames(cls):
        """
        Returns a list of XML output filenames that can be written by a calculation to the .save folder.

        Note that this includes all potential filenames across all known versions of Quantum ESPRESSO
        """
        return [cls._DATAFILE_XML_POST_6_2, cls._DATAFILE_XML_PRE_6_2]

    @abc.abstractmethod
    @classproperty
    def xml_filepaths(cls):
        """
        Returns a list of relative filepaths of XML files.

        This assumes that all xml filenames returned by the xml_filenames class property are written to
        the same .save output folder
        """
        pass

    # To be specified in the subclass:
    # _automatic_namelists = {
    #        'scf':   ['CONTROL', 'SYSTEM', 'ELECTRONS'],
    #        'nscf':  ['CONTROL', 'SYSTEM', 'ELECTRONS'],
    #        'bands': ['CONTROL', 'SYSTEM', 'ELECTRONS'],
    #        'relax': ['CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS'],
    #        'md':    ['CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS'],
    #        'vc-md':    ['CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS', 'CELL'],
    #        'vc-relax': ['CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS', 'CELL'],
    #        }

    # Keywords that cannot be set
    # If the length of the tuple is three, the third value is the value that
    # will be automatically set.
    # Note that some values (ibrav, nat, ntyp, ...) are overridden anyway
    #     _blocked_keywords = [('CONTROL', 'pseudo_dir'), # set later
    #          ('CONTROL', 'outdir'),  # set later
    #          ('CONTROL', 'prefix'),  # set later
    #          ('SYSTEM', 'ibrav'),  # set later
    #          ('SYSTEM', 'celldm'),
    #          ('SYSTEM', 'nat'),  # set later
    #          ('SYSTEM', 'ntyp'),  # set later
    #          ('SYSTEM', 'a'), ('SYSTEM', 'b'), ('SYSTEM', 'c'),
    #     ('SYSTEM', 'cosab'), ('SYSTEM', 'cosac'), ('SYSTEM', 'cosbc'),
    # ]

    # _use_kpoints = False

    @classmethod
    def define(cls, spec):
        super(BasePwCpInputGenerator, cls).define(spec)
        spec.input('code', valid_type=orm.Code, help='')
        spec.input('structure', valid_type=orm.StructureData, help='')
        spec.input('parameters', valid_type=orm.Dict, help='')
        spec.input('settings', valid_type=orm.Dict, required=False, help='')
        spec.input('parent_folder', valid_type=orm.RemoteData, required=False, help='')
        spec.input('vdw_table', valid_type=orm.SinglefileData, required=False, help='')
        spec.input_namespace('pseudos', valid_type=orm.UpfData, dynamic=True, help='')

    def prepare_for_submission(self, folder):
        """Create the input files from the input nodes passed to this instance of the `CalcJob`.

        :param folder: an `aiida.common.folders.Folder` to temporarily write files on disk
        :return: `aiida.common.datastructures.CalcInfo` instance
        """
        if 'settings' in self.inputs:
            settings = self.inputs.settings.get_dict()
            settings_dict = _uppercase_dict(settings, dict_name='settings')
        else:
            settings = {}

        # Check that a pseudo potential was specified for each kind present in the `StructureData`
        kinds = [kind.name for kind in self.inputs.structure.kinds]
        if set(kinds) != set(self.inputs.pseudos.keys()):
            raise exceptions.InputValidationError(
                'Mismatch between the defined pseudos and the list of kinds of the structure.\n'
                'Pseudos: {};\nKinds: {}'.format(', '.join(list(self.inputs.pseudos.keys())), ', '.join(list(kinds))))

        local_copy_list = []
        remote_copy_list = []
        remote_symlink_list = []

        # Create the subfolder that will contain the pseudopotentials
        folder.get_subfolder(self._PSEUDO_SUBFOLDER, create=True)
        # Create the subfolder for the output data (sometimes Quantum ESPRESSO codes crash if the folder does not exist)
        folder.get_subfolder(self._OUTPUT_SUBFOLDER, create=True)

        # If present, add also the Van der Waals table to the pseudo dir. Note that the name of the table is not checked
        # but should be the one expected by Quantum ESPRESSO.
        if 'vdw_table' in self.inputs:
            src_path = self.inputs.vdw_table.get_file_abs_path()
            dst_path = os.path.join(self._PSEUDO_SUBFOLDER, os.path.split(self.inputs.vdw_table.get_file_abs_path())[1])
            local_copy_list.append((src_path, dst_path))

        if 'hubbard_file' in self.inputs:
            src_path = self.inputs.hubbard_file.get_file_abs_path()
            dst_path = self.input_file_name_hubbard_file
            local_copy_list.append((src_path, dst_path))

        arguments = [
            self.inputs.parameters,
            settings_dict,
            self.inputs.pseudos,
            self.inputs.structure,
            self.inputs.kpoints
        ]
        input_filecontent, local_copy_pseudo_list = self._generate_PWCPinputdata(*arguments)
        local_copy_list += local_copy_pseudo_list

        input_filename = folder.get_abs_path(self._INPUT_FILE_NAME)
        with io.open(input_filename, 'w') as handle:
            handle.write(input_filecontent)
