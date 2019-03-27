# -*- coding: utf-8 -*-
from __future__ import absolute_import
from aiida.common import exceptions
from aiida.orm.nodes.data.remote import RemoteData
from aiida_quantumespresso.calculations.namelists import NamelistsCalculation
from aiida_quantumespresso.calculations.pw import PwCalculation


class PpCalculation(NamelistsCalculation):
    """
    Pp.x code of the Quantum ESPRESSO distribution, handles the
    post-processing of charge-densities, potentials, ...
    For more information, refer to http://www.quantum-espresso.org/
    """
    def _init_internal_params(self):
        super(PpCalculation, self)._init_internal_params()
        self._default_namelists = ['INPUTPP', 'PLOT']
        self._FILPLOT = "aiida.filplot"
        self._blocked_keywords = [
            ('INPUTPP', 'outdir', self._OUTPUT_SUBFOLDER),
            ('INPUTPP', 'prefix', self._PREFIX),
            ('INPUTPP', 'filplot', self._FILPLOT),
        ]
        self._default_parser = None
        self._internal_retrieve_list = [self._FILPLOT]

    def use_parent_calculation(self, calc):
        """Set the parent calculation."""
        if not isinstance(calc, PwCalculation):
            raise ValueError('Parent calculation must be a PwCalculation')

        try:
            remote_folder = calc.get_outgoing(node_class=RemoteData, link_label_filter='remote_folder').one().node
        except ValueError:
            raise exceptions.UniquenessError('Parent calculation does not have a remote folder output node.')

        self.use_parent_folder(remote_folder)
