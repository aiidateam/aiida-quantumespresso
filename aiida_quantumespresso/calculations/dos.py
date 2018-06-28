# -*- coding: utf-8 -*-
import os
from aiida.common.utils import classproperty
from aiida.orm.data.folder import FolderData
from aiida.orm.data.remote import RemoteData
from aiida.orm.data.array.kpoints import KpointsData
from aiida_quantumespresso.calculations.namelists import NamelistsCalculation
from aiida_quantumespresso.calculations.pw import PwCalculation
from aiida_quantumespresso.calculations import BasePwCpInputGenerator


class DosCalculation(NamelistsCalculation):
    """
    Plugin for the dos.x code of the Quantum ESPRESSO distribution. Handles
    density of states calculations, and stores the resulting dos arrays and
    integrated dos arrays.
    For more information regarding dos.x
    refer to http://www.quantum-espresso.org/
    """

    def _init_internal_params(self):
        super(DosCalculation, self)._init_internal_params()

        self._DOS_FILENAME = 'aiida.dos'
        self._default_namelists = ['DOS']
        self._blocked_keywords = [
            ('DOS', 'fildos', self._DOS_FILENAME),
            ('DOS', 'outdir', self._OUTPUT_SUBFOLDER),
            ('DOS', 'prefix', self._PREFIX),
        ]
        self._internal_retrieve_list = [self._DOS_FILENAME]
        self._default_parser = 'quantumespresso.dos'

    def use_parent_calculation(self, calc):
        """
        Set the parent calculation,
        from which it will inherit the outputsubfolder.
        The link will be created from parent RemoteData and NamelistCalculation
        """
        if not isinstance(calc, PwCalculation):
            raise ValueError("Parent calculation must be a PwCalculation")
        from aiida.common.exceptions import UniquenessError
        localdatas = [_[1] for _ in calc.get_outputs(also_labels=True)]
        if len(localdatas) == 0:
            raise UniquenessError("No output retrieved data found in the parent"
                                  "calc, probably it did not finish yet, "
                                  "or it crashed")

        localdata = [_[1] for _ in calc.get_outputs(also_labels=True) if _[0] == 'remote_folder']
        localdata = localdata[0]
        self.use_parent_folder(localdata)
