# -*- coding: utf-8 -*-
from aiida.common import exceptions
from aiida.orm.data.remote import RemoteData
from aiida_quantumespresso.calculations.namelists import NamelistsCalculation
from aiida_quantumespresso.calculations.pw import PwCalculation


class ProjwfcCalculation(NamelistsCalculation):
    """
    Projwfc.x code of the Quantum ESPRESSO distribution, handles the the
    computation of projections of bloch wavefunctions onto atomic orbitals
    <Psi(n,k) | Y(theta,phi)R(r) >.
    For more information, refer to http://www.quantum-espresso.org/
    """
    def _init_internal_params(self):
        super(ProjwfcCalculation, self)._init_internal_params()
        # self._PROJWFC_FILENAME = 'aiida.pdos'
        self._default_namelists = ['PROJWFC']
        self._blocked_keywords = [
                                  ('PROJWFC','outdir',self._OUTPUT_SUBFOLDER),
                                  ('PROJWFC','prefix',self._PREFIX),
                                  ('PROJWFC','lsym',True),
                                  ('PROJWFC','lwrite_overlaps',False),
                                  ('PROJWFC','lbinary_data',False),
                                  ('PROJWFC','kresolveddos',False),
                                  ('PROJWFC','tdosinboxes',False),
                                  ('PROJWFC','plotboxes',False),

                                 ]
        self._default_parser = 'quantumespresso.projwfc'
        self._internal_retrieve_list = [self._PREFIX+".pdos*"]

    def use_parent_calculation(self, calc):
        """Set the parent calculation."""
        if not isinstance(calc, PwCalculation):
            raise ValueError('Parent calculation must be a PwCalculation')

        try:
            remote_folder = calc.get_outgoing(node_class=RemoteData, link_label_filter='remote_folder').one().node
        except ValueError:
            raise exceptions.UniquenessError('Parent calculation does not have a remote folder output node.')

        self.use_parent_folder(remote_folder)
