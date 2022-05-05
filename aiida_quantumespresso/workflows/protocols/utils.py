# -*- coding: utf-8 -*-
"""Utilities to manipulate the workflow input protocols."""
import pathlib
from typing import Optional, Union

from aiida.plugins import DataFactory, GroupFactory
import yaml

StructureData = DataFactory('core.structure')
PseudoPotentialFamily = GroupFactory('pseudo.family')


class ProtocolMixin:
    """Utility class for processes to build input mappings for a given protocol based on a YAML configuration file."""

    @classmethod
    def get_protocol_filepath(cls) -> pathlib.Path:
        """Return the ``pathlib.Path`` to the ``.yaml`` file that defines the protocols."""
        raise NotImplementedError

    @classmethod
    def get_default_protocol(cls) -> str:
        """Return the default protocol for a given workflow class.

        :param cls: the workflow class.
        :return: the default protocol.
        """
        return cls._load_protocol_file()['default_protocol']

    @classmethod
    def get_available_protocols(cls) -> dict:
        """Return the available protocols for a given workflow class.

        :param cls: the workflow class.
        :return: dictionary of available protocols, where each key is a protocol and value is another dictionary that
            contains at least the key `description` and optionally other keys with supplementary information.
        """
        data = cls._load_protocol_file()
        return {protocol: {'description': values['description']} for protocol, values in data['protocols'].items()}

    @classmethod
    def get_protocol_inputs(
        cls,
        protocol: Optional[dict] = None,
        overrides: Union[dict, pathlib.Path, None] = None,
    ) -> dict:
        """Return the inputs for the given workflow class and protocol.

        :param cls: the workflow class.
        :param protocol: optional specific protocol, if not specified, the default will be used
        :param overrides: dictionary of inputs that should override those specified by the protocol. The mapping should
            maintain the exact same nesting structure as the input port namespace of the corresponding workflow class.
        :return: mapping of inputs to be used for the workflow class.
        """
        data = cls._load_protocol_file()
        protocol = protocol or data['default_protocol']

        try:
            protocol_inputs = data['protocols'][protocol]
        except KeyError as exception:
            raise ValueError(
                f'`{protocol}` is not a valid protocol. Call ``get_available_protocols`` to show available protocols.'
            ) from exception
        inputs = recursive_merge(data['default_inputs'], protocol_inputs)
        inputs.pop('description')

        if isinstance(overrides, pathlib.Path):
            with overrides.open() as file:
                overrides = yaml.safe_load(file)

        if overrides:
            return recursive_merge(inputs, overrides)

        return inputs

    @classmethod
    def _load_protocol_file(cls) -> dict:
        """Return the contents of the protocol file for workflow class."""
        with cls.get_protocol_filepath().open() as file:
            return yaml.safe_load(file)


def recursive_merge(left: dict, right: dict) -> dict:
    """Recursively merge two dictionaries into a single dictionary.

    If any key is present in both ``left`` and ``right`` dictionaries, the value from the ``right`` dictionary is
    assigned to the key.

    :param left: first dictionary
    :param right: second dictionary
    :return: the recursively merged dictionary
    """
    import collections

    # Note that a deepcopy is not necessary, since this function is called recusively.
    right = right.copy()

    for key, value in left.items():
        if key in right:
            if isinstance(value, collections.abc.Mapping) and isinstance(right[key], collections.abc.Mapping):
                right[key] = recursive_merge(value, right[key])

    merged = left.copy()
    merged.update(right)

    return merged


def get_magnetization_parameters() -> dict:
    """Return the mapping of suggested initial magnetic moments for each element.

    :returns: the magnetization parameters.
    """
    with (pathlib.Path(__file__).resolve().parent / 'magnetization.yaml').open() as handle:
        return yaml.safe_load(handle)


def get_starting_magnetization(
    structure: StructureData,
    pseudo_family: PseudoPotentialFamily,
    initial_magnetic_moments: Optional[dict] = None
) -> dict:
    """Return the dictionary with starting magnetization for each kind in the structure.

    :param structure: the structure.
    :param pseudo_family: pseudopotential family.
    :param initial_magnetic_moments: dictionary mapping each kind in the structure to its magnetic moment.
    :returns: dictionary of starting magnetizations.
    """
    if initial_magnetic_moments is not None:

        nkinds = len(structure.kinds)

        if sorted(initial_magnetic_moments.keys()) != sorted(structure.get_kind_names()):
            raise ValueError(f'`initial_magnetic_moments` needs one value for each of the {nkinds} kinds.')

        return {
            kind.name: initial_magnetic_moments[kind.name] / pseudo_family.get_pseudo(element=kind.symbol).z_valence
            for kind in structure.kinds
        }

    magnetic_parameters = get_magnetization_parameters()
    starting_magnetization = {}

    for kind in structure.kinds:
        magnetic_moment = magnetic_parameters[kind.symbol]['magmom']

        if magnetic_moment == 0:
            magnetization = magnetic_parameters['default_magnetization']
        else:
            z_valence = pseudo_family.get_pseudo(element=kind.symbol).z_valence
            magnetization = magnetic_moment / float(z_valence)

        starting_magnetization[kind.name] = magnetization

    return starting_magnetization
