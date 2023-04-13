# -*- coding: utf-8 -*-
"""Utility functions for validation of command line interface parameter inputs."""
from aiida.cmdline.utils import decorators
import click


@decorators.with_dbenv()
def validate_kpoints_mesh(ctx, param, value):
    """Command line option validator for a kpoints mesh tuple.

    The value should be a tuple of three positive integers out of which a KpointsData object will be created with a mesh
    equal to the tuple.

    :param ctx: internal context of the click.command
    :param param: the click Parameter, i.e. either the Option or Argument to which the validator is hooked up
    :param value: a tuple of three positive integers
    :returns: a KpointsData instance
    """
    # pylint: disable=unused-argument
    from aiida.orm import KpointsData

    if not value:
        return None

    if any(not isinstance(integer, int) for integer in value) or any(int(i) <= 0 for i in value):
        raise click.BadParameter('all values of the tuple should be positive greater than zero integers')

    try:
        kpoints = KpointsData()
        kpoints.set_kpoints_mesh(value)
    except ValueError as exception:
        raise click.BadParameter(f'failed to create a KpointsData mesh out of {value}\n{exception}')

    return kpoints


@decorators.with_dbenv()
def validate_hubbard_parameters(structure, parameters, hubbard_u=None, hubbard_v=None, hubbard_file=None):
    """Validate Hubbard input parameters and update the parameters input node accordingly.

    :param structure: the StructureData node that will be used in the inputs
    :param parameters: the Dict node that will be used in the inputs
    :param hubbard_u: the Hubbard U inputs values from the cli
    :param hubbard_v: the Hubbard V inputs values from the cli
    :param hubbard_file: a SinglefileData with Hubbard parameters
    :raises ValueError: if the input is invalid
    """
    if len([value for value in [hubbard_u, hubbard_v, hubbard_file] if value]) > 1:
        raise ValueError('the hubbard_u, hubbard_v and hubbard_file options are mutually exclusive')

    if hubbard_file:

        parameters['SYSTEM']['lda_plus_u'] = True
        parameters['SYSTEM']['lda_plus_u_kind'] = 2
        parameters['SYSTEM']['hubbard_parameters'] = 'file'

    elif hubbard_v:

        parameters['SYSTEM']['lda_plus_u'] = True
        parameters['SYSTEM']['lda_plus_u_kind'] = 2
        parameters['SYSTEM']['hubbard_parameters'] = 'input'
        parameters['SYSTEM']['hubbard_v'] = []

        for value in hubbard_v:
            parameters['SYSTEM']['hubbard_v'].append(value)

    elif hubbard_u:

        structure_kinds = structure.get_kind_names()
        hubbard_kinds = [value[0] for value in hubbard_u]

        if not set(hubbard_kinds).issubset(structure_kinds):
            raise ValueError('kinds in the specified Hubbard U is not a strict subset of the structure kinds')

        parameters['SYSTEM']['lda_plus_u'] = True
        parameters['SYSTEM']['lda_plus_u_kind'] = 0
        parameters['SYSTEM']['hubbard_u'] = {}

        for kind, value in hubbard_u:
            parameters['SYSTEM']['hubbard_u'][kind] = value


def validate_starting_magnetization(structure, parameters, starting_magnetization=None):
    """Validate starting magnetization parameters and update the parameters input node accordingly.

    :param structure: the StructureData node that will be used in the inputs
    :param parameters: the Dict node that will be used in the inputs
    :param starting_magnetization: the starting magnetization inputs values from the cli
    :raises ValueError: if the input is invalid
    """
    if not starting_magnetization:
        return

    structure_kinds = structure.get_kind_names()
    magnetization_kinds = [kind for kind, magnetization in starting_magnetization]

    if not set(magnetization_kinds).issubset(structure_kinds):
        raise ValueError('kinds in the specified starting magnetization is not a strict subset of the structure kinds')

    parameters['SYSTEM']['nspin'] = 2
    parameters['SYSTEM']['starting_magnetization'] = {}

    for kind, magnetization in starting_magnetization:
        parameters['SYSTEM']['starting_magnetization'][kind] = magnetization


def validate_smearing(parameters, smearing=None):
    """Validate smearing parameters and update the parameters input node accordingly.

    :param parameters: the Dict node that will be used in the inputs
    :param smearing: a tuple of a string and float corresponding to type of smearing and the degauss value
    :raises ValueError: if the input is invalid
    """
    if not any(smearing):
        return

    valid_smearing_types = {
        'gaussian': ['gaussian', 'gauss'],
        'methfessel-paxton': ['methfessel-paxton', 'm-p', 'mp'],
        'marzari-vanderbilt': ['marzari-vanderbilt', 'cold', 'm-v', 'mv'],
        'fermi-dirac': ['fermi-dirac', 'f-d', 'fd'],
    }

    for _, options in valid_smearing_types.items():
        if smearing[0] in options:
            break
    else:
        raise ValueError(
            f'the smearing type "{smearing[0]}" is invalid, choose from {", ".join(list(valid_smearing_types.keys()))}'
        )

    if not isinstance(smearing[1], float):
        raise ValueError('the smearing value should be a float')

    parameters['SYSTEM']['occupations'] = 'smearing'
    parameters['SYSTEM']['smearing'] = smearing[0]
    parameters['SYSTEM']['degauss'] = smearing[1]
