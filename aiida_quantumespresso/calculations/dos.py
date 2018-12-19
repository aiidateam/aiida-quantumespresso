# -*- coding: utf-8 -*-
from aiida.common import exceptions
from aiida.orm.data.remote import RemoteData
from aiida_quantumespresso.calculations.namelists import NamelistsCalculation
from aiida_quantumespresso.calculations.pw import PwCalculation


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
        """Set the parent calculation."""
        if not isinstance(calc, PwCalculation):
            raise ValueError('Parent calculation must be a PwCalculation')

        try:
            remote_folder = calc.get_outgoing(node_class=RemoteData, link_label_filter='remote_folder').one().node
        except ValueError:
            raise exceptions.UniquenessError('Parent calculation does not have a remote folder output node.')

        self.use_parent_folder(remote_folder)
