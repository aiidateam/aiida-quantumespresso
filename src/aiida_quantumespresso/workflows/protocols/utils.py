"""Utilities to manipulate the workflow input protocols."""

import copy
import pathlib
import warnings
from typing import Optional, Union

import yaml
from aiida.common.warnings import AiidaDeprecationWarning
from aiida.orm import StructureData
from aiida_pseudo.groups.family import PseudoPotentialFamily
from plumpy import PortNamespace

from aiida_quantumespresso.common.types import SpinType


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
            alias_protocol = cls._check_if_alias(protocol)
            if alias_protocol is not None:
                protocol_inputs = data['protocols'][alias_protocol]
            else:
                raise ValueError(
                    f'`{protocol}` is not a valid protocol. Call ``get_available_protocols`` to show available '
                    'protocols.'
                ) from exception
        inputs = recursive_merge(data['default_inputs'], protocol_inputs)
        inputs.pop('description')

        if isinstance(overrides, pathlib.Path):
            with overrides.open() as file:
                overrides = yaml.safe_load(file)

        if overrides:
            cls._validate_override_keys(overrides)
            return recursive_merge(inputs, overrides)

        return inputs

    @staticmethod
    def set_default_resources(options: dict, scheduler_type: str) -> dict:
        """Set default resources based on the scheduler type of a computer.

        Temporary workaround to keep the default resources for the direct and Slurm schedulers until defaults can be
        properly configured on computers, see https://github.com/aiidateam/aiida-quantumespresso/pull/1011

        :param options: the options dictionary.
        :param scheduler_type: the scheduler type, e.g. `core.slurm`.
        """
        new_options = copy.deepcopy(options)

        if scheduler_type in (
            'core.direct',
            'core.slurm',
            'core.pbspro',
            'core.torque',
        ):
            new_options.setdefault('resources', {}).setdefault('num_machines', 1)
        if scheduler_type in ('core.sge',):
            new_options.setdefault('resources', {}).setdefault('parallel_env', 'mpi')
            new_options.setdefault('resources', {}).setdefault('tot_num_mpiprocs', 1)

        return new_options

    @classmethod
    def _load_protocol_file(cls) -> dict:
        """Return the contents of the protocol file for workflow class."""
        with cls.get_protocol_filepath().open() as file:
            return yaml.safe_load(file)

    @staticmethod
    def _check_if_alias(alias: str):
        """Check if a given alias corresponds to a valid protocol."""
        aliases_dict = {
            'moderate': 'balanced',
            'precise': 'stringent',
        }
        return aliases_dict.get(alias)

    @classmethod
    def _validate_override_keys(cls, overrides):
        """Validate that override keys match either the process spec or protocol inputs.

        Recursively checks the `overrides` dictionary against the input port namespace of the work chain, extended
        with protocol inputs, and emits warnings for any unrecognised keys.
        """

        def port_namespace_to_dict(namespace: PortNamespace):
            """Recursively convert a PortNamespace into a nested dict structure."""
            return {
                key: port_namespace_to_dict(port) if isinstance(port, PortNamespace) else None
                for key, port in namespace.items()
            }

        def recursive_key_check(inputs_mapping: dict, overrides: dict, path=''):
            """Recursively check that all the provided keys in the `overrides` are in the `inputs_mapping`."""

            for key, value in overrides.items():
                full_key = f'{path}.{key}' if path else key
                if key not in inputs_mapping:
                    warnings.warn(f'Found unrecognised key in overrides: {full_key}')
                    continue
                if isinstance(value, dict) and isinstance(inputs_mapping.get(key), dict):
                    recursive_key_check(inputs_mapping[key], value, full_key)

        inputs_mapping = recursive_merge(cls.get_protocol_inputs(), port_namespace_to_dict(cls.spec().inputs))
        recursive_key_check(inputs_mapping, overrides)


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
        if (
            key in right
            and isinstance(value, collections.abc.Mapping)
            and isinstance(right[key], collections.abc.Mapping)
        ):
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


def get_magnetization(
    structure: StructureData,
    z_valences: dict,
    initial_magnetic_moments: Optional[dict] = None,
    spin_type: SpinType = SpinType.COLLINEAR,
) -> dict:
    """Return a magnetization dictionary with for each kind in the structure.

    The returned dictionary always has three keys, corresponding to the Quantum ESPRESSO inputs:

    * `starting_magnetization`
    * `angle1`
    * `angle2`

    In case the `spin_type` is set to `SpinType.COLLINEAR`, the values for `angle1` and `angle2` will be set to `None`.

    :param structure: the structure.
    :param z_valences: dictionary mapping each kind in the structure to the number of valence electrons in the pseudo
        potential.
    :param initial_magnetic_moments: dictionary mapping each kind in the structure to its magnetic moment.
    :param spin_type: the `SpinType` of the calculation.
    :returns: dictionary of the magnetization.
    """
    magnetization = {
        'starting_magnetization': {},
        'angle1': {} if spin_type in [SpinType.NON_COLLINEAR, SpinType.SPIN_ORBIT] else None,
        'angle2': {} if spin_type in [SpinType.NON_COLLINEAR, SpinType.SPIN_ORBIT] else None,
    }

    if sorted(z_valences.keys()) != sorted(structure.get_kind_names()):
        raise ValueError(f'`z_valences` needs one value for each of the {len(structure.kinds)} kinds.')

    if initial_magnetic_moments is not None:
        if sorted(initial_magnetic_moments.keys()) != sorted(structure.get_kind_names()):
            raise ValueError(
                f'`initial_magnetic_moments` needs one value for each of the {len(structure.kinds)} kinds.'
            )

        for kind in structure.kinds:
            magmom = initial_magnetic_moments[kind.name]

            if isinstance(magmom, (int, float)):
                scaled_magmom = magmom / z_valences[kind.name]
                magnetization['starting_magnetization'][kind.name] = scaled_magmom

                if spin_type in [SpinType.NON_COLLINEAR, SpinType.SPIN_ORBIT]:
                    magnetization['angle1'][kind.name] = 0.0
                    magnetization['angle2'][kind.name] = 0.0

            elif isinstance(magmom, (list, tuple)):
                if spin_type not in [SpinType.NON_COLLINEAR, SpinType.SPIN_ORBIT]:
                    raise TypeError(
                        f'Spin type is set to `{spin_type}` but a `{type(magmom)}` is provided for the magnetic '
                        f'moment of kind `{kind.name}`.'
                    )

                scaled_magmom = magmom[0] / z_valences[kind.name]
                magnetization['starting_magnetization'][kind.name] = scaled_magmom
                magnetization['angle1'][kind.name] = magmom[1]
                magnetization['angle2'][kind.name] = magmom[2]
            else:
                raise TypeError(f'Unrecognised type for magnetic moment of kind `{kind.name}`: {type(magmom)}.')

        return magnetization

    # The following block deals with the case where the input `structure` is a `MagneticStructureData`, which is
    # implemented in aiida_wannier90_workflows. Here we can read the magnetic moment from the structure.
    # TODO: Deprecate/Remove this block when we provide support for the new `StructureData` instead.
    if hasattr(structure, 'has_magmom'):
        collinear = structure.is_collin_mag()

        if not collinear and spin_type == SpinType.COLLINEAR:
            raise ValueError(
                'Input `MagneticStructureData` has non-collinear magnetism defined but `spin_type` is set to COLLINEAR.'
            )

        for kind in structure.kinds:
            magmom = kind.get_magmom_coord()

            scaled_magmom = magmom[0] / z_valences[kind.name]
            magnetization['starting_magnetization'][kind.name] = scaled_magmom

            if spin_type in [SpinType.NON_COLLINEAR, SpinType.SPIN_ORBIT]:
                magnetization['angle1'][kind.name] = magmom[1]
                magnetization['angle2'][kind.name] = magmom[2]

        return magnetization

    # End of `MagneticStructureData` block

    magnetic_parameters = get_magnetization_parameters()

    for kind in structure.kinds:
        magnetic_moment = magnetic_parameters[kind.symbol]['magmom']

        magnetization['starting_magnetization'][kind.name] = (
            magnetic_parameters['default_magnetization']
            if magnetic_moment == 0
            else magnetic_moment / z_valences[kind.name]
        )
        if spin_type in [SpinType.NON_COLLINEAR, SpinType.SPIN_ORBIT]:
            magnetization['angle1'][kind.name] = 0.0
            magnetization['angle2'][kind.name] = 0.0

    return magnetization


def get_starting_magnetization(
    structure: StructureData,
    pseudo_family: PseudoPotentialFamily,
    initial_magnetic_moments: Optional[dict] = None,
) -> dict:
    """Return the dictionary with starting magnetization for each kind in the structure.

    :param structure: the structure.
    :param pseudo_family: pseudopotential family.
    :param initial_magnetic_moments: dictionary mapping each kind in the structure to its magnetic moment.
    :returns: dictionary of starting magnetizations.
    """
    warnings.warn(
        '`get_starting_magnetization` is deprecated, use `get_magnetization` instead.',
        AiidaDeprecationWarning,
    )

    if initial_magnetic_moments is not None:
        nkinds = len(structure.kinds)

        if sorted(initial_magnetic_moments.keys()) != sorted(structure.get_kind_names()):
            raise ValueError(f'`initial_magnetic_moments` needs one value for each of the {nkinds} kinds.')

        return {
            kind.name: initial_magnetic_moments[kind.name] / pseudo_family.get_pseudo(element=kind.symbol).z_valence
            for kind in structure.kinds
        }

    starting_magnetization = {}
    try:
        structure.has_magmom()
    except AttributeError:
        # Normal StructureData, no magmom in structure
        magnetic_parameters = get_magnetization_parameters()

        for kind in structure.kinds:
            magnetic_moment = magnetic_parameters[kind.symbol]['magmom']

            if magnetic_moment == 0:
                magnetization = magnetic_parameters['default_magnetization']
            else:
                z_valence = pseudo_family.get_pseudo(element=kind.symbol).z_valence
                magnetization = magnetic_moment / float(z_valence)

            starting_magnetization[kind.name] = magnetization
    else:
        # MagneticStructureData, currently implemented in aiida_wannier90_workflows.
        # Read magmom from structure<MagneticStructureData>
        collinear = structure.is_collin_mag()
        for kind in structure.kinds:
            magmom = kind.get_magmom_coord()[0] if collinear else kind.get_magmom_coord(coord='cartesian')[2]
            starting_magnetization[kind.name] = magmom / pseudo_family.get_pseudo(element=kind.symbol).z_valence

    return starting_magnetization


def get_starting_magnetization_noncolin(
    structure: StructureData,
    pseudo_family: PseudoPotentialFamily,
    initial_magnetic_moments: Optional[dict] = None
) -> tuple:
    """Return the dictionary with starting magnetization for each kind in the structure.

    :param structure: the structure.
    :param pseudo_family: pseudopotential family.
    :param initial_magnetic_moments: dictionary mapping each kind in the structure to its magnetic moment.
    :returns: dictionary of starting magnetizations.
    """
    # try:
    #     structure.mykinds
    # except AttributeError:
    #     raise TypeError(f"structure<{structure.pk}> do not have magmom")
    starting_magnetization = {}
    angle1 = {}
    angle2 = {}

    if initial_magnetic_moments is not None:

        nkinds = len(structure.kinds)

        if sorted(initial_magnetic_moments.keys()) != sorted(structure.get_kind_names()):
            raise ValueError(f'`initial_magnetic_moments` needs one value for each of the {nkinds} kinds.')

        for kind in structure.kinds:
            magmom = initial_magnetic_moments[kind.name]
            if isinstance(magmom, Union[int, float]):
                starting_magnetization[kind.name] = magmom / pseudo_family.get_pseudo(element=kind.symbol).z_valence
                angle1[kind.name] = 0.0
                angle2[kind.name] = 0.0
            else:  # tuple of 3 float (r, theta, phi)
                starting_magnetization[kind.name
                                       ] = 2 * magmom[0] / pseudo_family.get_pseudo(element=kind.symbol).z_valence
                angle1[kind.name] = magmom[1]
                angle2[kind.name] = magmom[2]
    try:
        structure.mykinds
    except AttributeError:
        # Normal StructureData, no magmom in structure
        magnetic_parameters = get_magnetization_parameters()

        for kind in structure.kinds:
            magnetic_moment = magnetic_parameters[kind.symbol]['magmom']

            if magnetic_moment == 0:
                magnetization = magnetic_parameters['default_magnetization']
            else:
                z_valence = pseudo_family.get_pseudo(element=kind.symbol).z_valence
                magnetization = magnetic_moment / float(z_valence)

            starting_magnetization[kind.name] = magnetization
            angle1[kind.name] = 0.0
            angle2[kind.name] = 0.0
    else:
        # Self defined myStructureData, read magmom from structure
        for kind in structure.mykinds:
            magmom = kind.get_magmom_coord()
            starting_magnetization[kind.name] = 2 * magmom[0] / pseudo_family.get_pseudo(element=kind.symbol).z_valence
            angle1[kind.name] = magmom[1]
            angle2[kind.name] = magmom[2]

    return starting_magnetization, angle1, angle2
