# -*- coding: utf-8 -*-
"""Utility class for handling the :class:`aiida_quantumespresso.data.hubbard_structure.HubbardStructureData`."""
# pylint: disable=no-name-in-module

from aiida import orm
from aiida.common.exceptions import MissingEntryPointError
from aiida.engine import calcfunction
from aiida.orm import StructureData as LegacyStructureData
from aiida.plugins import DataFactory

try:
    StructureData = DataFactory('atomistic.structure')
except MissingEntryPointError:
    structures_classes = (LegacyStructureData,)
else:
    structures_classes = (LegacyStructureData, StructureData)


class MagneticUtils:  # pylint: disable=too-few-public-methods
    """Class to manage the magnetic structure of the atomistic `LegacyStructureData`.

    It contains methods to manipulate the magne tic structure in such a way to produce
    the correct input for QuantumESPRESSO calculations.
    """

    def __init__(
        self,
        structure: structures_classes,
    ):
        """Set a the `StructureData` to manipulate."""
        if isinstance(structure, StructureData):
            if 'magmoms' not in structure.get_defined_properties():
                raise ValueError('The input structure does not contain magnetic moments.')
            self.structure = structure
        else:
            raise ValueError('input is not of type atomistic `StructureData')

    def generate_magnetic_namelist(self, parameters):
        """Generate the magnetic namelist for Quantum ESPRESSO.

        :param parameters: dictionary of inputs for the Quantum ESPRESSO calculation.
        """
        if 'nspin' not in parameters['SYSTEM'] and 'noncolin' not in parameters['SYSTEM']:
            raise ValueError("The input parameters must contain the 'nspin' or the 'noncolin' key.")

        namelist = {'starting_magnetization': {}, 'angle1': {}, 'angle2': {}}

        if parameters['SYSTEM'].get('nspin', None) == 2:
            namelist.pop('angle1')
            namelist.pop('angle2')
            if self.structure.is_collinear:
                for kind, magmom in zip(self.structure.kinds, self.structure.magmoms):
                    # this should be fixed, now only magmom_z is considered...
                    if magmom[2] != 0:
                        namelist['starting_magnetization'][kind] = magmom[2]
            else:
                raise NotImplementedError(
                    'The input structure is not collinear, but you choose collinear calculations.'
                )
        elif parameters['SYSTEM']['noncolin']:
            for site in self.structure.sites:
                for variable, value in namelist.items():
                    value[site.kinds] = site.get_magmom_coord(coord='spherical')[variable]

        return namelist


@calcfunction
def generate_structure_with_magmoms(input_structure: structures_classes, input_magnetic_moments: orm.List):
    """Generate a new structure with the magnetic moments for each site.

    :param input_structure: the input structure to add the magnetic moments.
    :param input_magnetic_moments: the magnetic moments for each site, represented as a float (see below).

    For now, only supports collinear magnetic moments, i.e. atomic magnetizations (along z axis).
    """
    magmoms = input_magnetic_moments.get_list()
    if len(input_structure.sites) != len(magmoms):
        raise ValueError('The input structure and the magnetic moments must have the same length.')

    mutable_structure = input_structure.get_value()
    mutable_structure.clear_sites()
    for site, magmom in zip(input_structure.sites, magmoms):
        mutable_structure.append_atom(
            **{
                'position': site.position,
                'symbol': site.symbol,
                'kind': site.kind,
                'weight': site.weight,
                'magmom': [0, 0, magmom] if isinstance(magmom, float) else magmom  # 3D vector
            }
        )

    output_structure = StructureData.from_mutable(mutable_structure, detect_kinds=True)

    return output_structure
