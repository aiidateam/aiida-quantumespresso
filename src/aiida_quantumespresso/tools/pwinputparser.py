# -*- coding: utf-8 -*-
"""Utilities to parse Quantum ESPRESSO pw.x input files into AiiDA nodes or builders."""
import copy
import re

from aiida.common.folders import Folder
from aiida.orm import Dict, load_code
from aiida.plugins import CalculationFactory, DataFactory
import numpy as np
from qe_tools.parsers import PwInputFile as BasePwInputFile

from .base import StructureParseMixin

UpfData = DataFactory('pseudo.upf')


class PwInputFile(StructureParseMixin, BasePwInputFile):
    """Parser of Quantum ESPRESSO pw.x input file into AiiDA nodes.

    .. note:: This mixes in :class:`~aiida_quantumespresso.tools.base.StructureParseMixin` which adds the functionality
        to parse a :class:`~aiida.orm.nodes.data.structure.StructureData` from the input file, instead of a plain
        dictionary returned by ``qe_tools.parsers.qeinputparser.get_structure_from_qeinput``. Note that one cannot
        directly add this functionality to a sub class of ``~qe_tools.parsers.qeinputparser.QeInputFile`` and then
        subsequently sub class that here, because the ``~qe_tools.parsers.qeinputparser.CpInputFile`` is also
        required and sub classing both leads to problems with the MRO.
    """

    def get_kpointsdata(self):
        """Return a `KpointsData` object based on the data in the input file.

        .. note:: If the calculation uses only the gamma k-point (`if self.k_points['type'] == 'gamma'`), it is
            necessary to also attach a settings node to the calculation with `gamma_only = True`.

        :return: KpointsData object of the kpoints in the input file
        :rtype: :class:`~aiida.orm.nodes.data.array.kpoints.KpointsData`
        :raises NotImplementedError: if the kpoints are in a format not yet supported.
        """
        from aiida.orm.nodes.data.array.kpoints import KpointsData

        kpoints = KpointsData()
        structure = self.get_structuredata()
        kpoints.set_cell_from_structure(structure)

        # Set the kpoints and weights, doing any necessary units conversion.
        if self.k_points['type'] == 'crystal':  # relative to reciprocal lattice vectors
            kpoints.set_kpoints(self.k_points['points'], weights=self.k_points['weights'])
        elif self.k_points['type'] == 'tpiba':  # Cartesian; units of 2*pi/alat
            alat = np.linalg.norm(structure.cell[0])  # alat in Angstrom
            kpoints.set_kpoints(
                np.array(self.k_points['points']) * (2. * np.pi / alat),
                weights=self.k_points['weights'],
                cartesian=True
            )
        elif self.k_points['type'] == 'automatic':
            kpoints.set_kpoints_mesh(self.k_points['points'], offset=self.k_points['offset'])
        elif self.k_points['type'] == 'gamma':
            kpoints.set_kpoints_mesh([1, 1, 1])
        else:
            raise NotImplementedError(f"support for units {self.k_points['type']} not yet implemented")

        return kpoints


def create_builder_from_file(input_folder, input_file_name, code, metadata, pseudo_folder_path=None):
    """Create a populated process builder for a `PwCalculation` from a standard QE input file and pseudo (upf) files.

    :param input_folder: the folder containing the input file
    :type input_folder: aiida.common.folders.Folder or str
    :param input_file_name: the name of the input file
    :type input_file_name: str
    :param code: the code associated with the calculation
    :type code: aiida.orm.AbstractCode or str
    :param metadata: metadata values for the calculation (e.g. resources)
    :type metadata: dict
    :param pseudo_folder_path: the folder containing the upf files (if None, then input_folder is used)
    :type pseudo_folder_path: aiida.common.folders.Folder or str or None
    :raises NotImplementedError: if the structure is not ibrav=0
    :return: a builder instance for PwCalculation
    """
    PwCalculation = CalculationFactory('quantumespresso.pw')

    builder = PwCalculation.get_builder()
    builder.metadata = metadata

    if isinstance(code, str):
        code = load_code(code)
    builder.code = code

    # read input_file
    if isinstance(input_folder, str):
        input_folder = Folder(input_folder)

    with input_folder.open(input_file_name) as input_file:
        parsed_file = PwInputFile(input_file.read())

    builder.structure = parsed_file.get_structuredata()
    builder.kpoints = parsed_file.get_kpointsdata()

    # Then, strip the namelist items that the plugin doesn't allow or sets later.
    # NOTE: If any of the position or cell units are in alat or crystal
    # units, that will be taken care of by the input parsing tools, and
    # we are safe to fake that they were never there in the first place.
    parameters_dict = copy.deepcopy(parsed_file.namelists)
    for namelist, blocked_key in PwCalculation._blocked_keywords:  # pylint: disable=protected-access
        for key in list(parameters_dict[namelist].keys()):
            # take into account that celldm and celldm(*) must be blocked
            if re.sub('[(0-9)]', '', key) == blocked_key:
                parameters_dict[namelist].pop(key, None)
    builder.parameters = Dict(parameters_dict)

    # Get or create a UpfData node for the pseudopotentials used for the calculation.
    pseudos_map = {}
    if pseudo_folder_path is None:
        pseudo_folder_path = input_folder
    if isinstance(pseudo_folder_path, str):
        pseudo_folder_path = Folder(pseudo_folder_path)
    names = parsed_file.atomic_species['names']
    pseudo_file_names = parsed_file.atomic_species['pseudo_file_names']
    pseudo_file_map = {}
    for name, fname in zip(names, pseudo_file_names):
        if fname not in pseudo_file_map:
            local_path = pseudo_folder_path.get_abs_path(fname)
            with open(local_path, 'rb') as handle:
                upf = UpfData(handle)
            pseudo_file_map[fname] = upf
        pseudos_map[name] = pseudo_file_map[fname]
    builder.pseudos = pseudos_map

    settings_dict = {}
    if parsed_file.k_points['type'] == 'gamma':
        settings_dict['gamma_only'] = True

    # If there are any fixed coordinates (i.e. force modification) present in the input file, specify in settings
    fixed_coords = parsed_file.atomic_positions['fixed_coords']
    # Function ``any()`` only works for 1-dimensional lists so we have to call it twice manually.
    if any((any(fc_xyz) for fc_xyz in fixed_coords)):
        settings_dict['FIXED_COORDS'] = fixed_coords

    if settings_dict:
        builder.settings = settings_dict

    return builder
