# -*- coding: utf-8 -*-
from __future__ import absolute_import
import os
from aiida.common import exceptions
from aiida.common.lang import classproperty
from aiida.orm.nodes.data.remote import RemoteData
from aiida.orm.nodes.data.array.kpoints import KpointsData
from aiida_quantumespresso.calculations.namelists import NamelistsCalculation
from aiida_quantumespresso.calculations.q2r import Q2rCalculation


class MatdynCalculation(NamelistsCalculation):
    """
    matdyn.x code of the Quantum ESPRESSO distribution, used to obtain the
    phonon frequencies in reciprocal space from the interatomic force constants given by q2r.
    For more information, refer to http://www.quantum-espresso.org/
    """

    def _init_internal_params(self):
        super(MatdynCalculation, self)._init_internal_params()
                
        self._PHONON_FREQUENCIES_NAME = 'phonon_frequencies.dat'
        self._PHONON_MODES_NAME = 'phonon_displacements.dat'
        self._PHONON_DOS_NAME = 'phonon_dos.dat'
    
        self._default_namelists = ['INPUT']   
        
        self._blocked_keywords = [('INPUT','flfrq',self._PHONON_FREQUENCIES_NAME), # output freq.
                                  ('INPUT','flvec',self._PHONON_MODES_NAME), # output displ.
                                  ('INPUT','fldos',self._PHONON_DOS_NAME), # output dos
                                  ('INPUT','q_in_cryst_coord',True), # kpoints always in crystal coordinates
                                  # this is dynamically added in the _prepare_for_submission
                                  #('INPUT','flfrc',Q2rCalculation.FORCE_CONSTANTS_NAME), # input
                                 ]
    
        self._internal_retrieve_list = [self._PHONON_FREQUENCIES_NAME, 
                                        self._PHONON_DOS_NAME]
        
        # Default Matdyn output parser provided by AiiDA
        self._default_parser = 'quantumespresso.matdyn'
    
    @classproperty
    def _use_methods(cls):
        """
        Use_* methods for the matdyn class.
        """
        retdict = NamelistsCalculation._use_methods
        retdict.update({
            "kpoints": {
               'valid_types': (KpointsData),
               'additional_parameter': None,
               'linkname': 'kpoints',
               'docstring': ("Kpoints on which to calculate the phonon "
                             "frequencies"),
               },                        
            })
        return retdict

    def use_parent_calculation(self, calc):
        """Set the parent calculation."""
        if not isinstance(calc, Q2rCalculation):
            raise ValueError('Parent calculation must be a Q2rCalculation')

        try:
            remote_folder = calc.get_outgoing(node_class=RemoteData, link_label_filter='remote_folder').one().node
        except ValueError:
            raise exceptions.UniquenessError('Parent calculation does not have a remote folder output node.')

        self.use_parent_folder(remote_folder)
   
    def _get_following_text(self, inputdict, settings):
        """
        Add the kpoints after the namelist.
        
        This function should consume the content of inputdict (if it requires
        a different node) or the keys inside settings, using the 'pop' method,
        so that inputdict and settings should remain empty at the end of 
        _prepare_for_submission, if all flags/nodes were recognized
        """
        from aiida.common.exceptions import InputValidationError
        
        try:
            kpoints = inputdict.pop(self.get_linkname('kpoints'))
        except KeyError:
            raise InputValidationError("No kpoints specified for this calculation")
        if not isinstance(kpoints, KpointsData):
            raise InputValidationError("kpoints is not of type KpointsData")
        
        try: 
            klist = kpoints.get_kpoints()
        except AttributeError:
            klist = kpoints.get_kpoints_mesh(print_list=True)
        
        retlist = ["{}".format(len(klist))]
        for k in klist:
            retlist.append("{:18.10f} {:18.10f} {:18.10f}".format(*k))
        
        return "\n".join(retlist)+"\n"
    
    def _prepare_for_submission(self,tempfolder, inputdict): 
        from aiida.orm.nodes.data.singlefile import SinglefileData
        
        parent_calc_folder = inputdict.get(self.get_linkname('parent_folder'),
                                           None)
        
        if isinstance(parent_calc_folder, SinglefileData):
            self._blocked_keywords.append(
                ('INPUT', 'flfrc', os.path.split(
                    parent_calc_folder.get_file_abs_path())[1]  ))
        else:
            raise NotImplementedError(
                "Input different from SinglefileData is not supported"
                " yet for MatdynCalculation; it is {}".format(
                type(parent_calc_folder)))
            self._blocked_keywords.append(
                ('INPUT', 'flfrc',  Q2rCalculation._FORCE_CONSTANTS_NAME ))
                
        calcinfo = super(MatdynCalculation, self)._prepare_for_submission(
            tempfolder, inputdict)
        return calcinfo
        
