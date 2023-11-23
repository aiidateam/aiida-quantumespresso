# -*- coding: utf-8 -*-
"""Tools for nodes created by running the `PwCalculation` class."""
from aiida.common import AttributeDict, exceptions
from aiida.tools.calculations.base import CalculationTools

from aiida_quantumespresso.calculations.functions.create_magnetic_configuration import create_magnetic_configuration


class PwCalculationTools(CalculationTools):
    """Calculation tools for `PwCalculation`.

    Methods implemented here are available on any `CalcJobNode` produced by the `PwCalculation class through the `tools`
    attribute.
    """

    # pylint: disable=too-few-public-methods

    def get_scf_accuracy(self, index=0):
        """Return the array of SCF accuracy values for a given SCF cycle.

        :param index: the zero-based index of the desired SCF cycle
        :return: a list of SCF accuracy values of a certain SCF cycle.
        :raises ValueError: if the node does not have the `output_trajectory` output
        :raises ValueError: if `output_trajectory` does not have the `scf_accuracy` or `scf_iterations` arrays
        :raises IndexError: if the `index` is out of range
        """
        try:
            trajectory = self._node.outputs.output_trajectory
        except exceptions.NotExistent as exc:
            raise ValueError(f'{self._node} does not have the `output_trajectory` output node') from exc

        try:
            scf_accuracy = trajectory.get_array('scf_accuracy')
        except KeyError as exc:
            raise ValueError(f'{trajectory} does not contain the required `scf_accuracy` array') from exc

        try:
            scf_iterations = trajectory.get_array('scf_iterations')
        except KeyError as exc:
            raise ValueError(f'{trajectory} does not contain the required `scf_iterations` array') from exc

        number_of_frames = len(scf_iterations)

        if index < -number_of_frames or index >= number_of_frames:
            raise IndexError(f'invalid index {index}, must be between {-number_of_frames} and {number_of_frames - 1}')

        # building an scf_accuracy_index for easier manipulation
        scf_accuracy_index = [0]
        for i in scf_iterations:
            scf_accuracy_index.append(scf_accuracy_index[-1] + i)

        if index < 0:
            return scf_accuracy[scf_accuracy_index[index - 1]:scf_accuracy_index[index]]

        return scf_accuracy[scf_accuracy_index[index]:scf_accuracy_index[index + 1]]

    def get_magnetic_configuration(self, atol: float = 0.5, ztol: float = 0.05) -> AttributeDict:
        """Get the final magnetic configuration of a ``pw.x`` calculation."""
        try:
            structure = self._node.outputs.output_structure
        except AttributeError:
            structure = self._node.inputs.structure

        try:
            magnetic_moments = self._node.outputs.output_trajectory.get_array('atomic_magnetic_moments')[-1].tolist()
        except KeyError as exc:
            raise KeyError('Could not find the `atomic_magnetic_moments` in the `output_trajectory`.') from exc

        results = create_magnetic_configuration(
            structure,
            magnetic_moments,
            atol=atol,
            ztol=ztol,
            metadata={'store_provenance': False},
        )
        structure_kindname_position = [(site.kind_name, site.position) for site in structure.sites]
        allo_kindname_position = [(site.kind_name, site.position) for site in results['structure'].sites]
        structure_kindnames_sorted = sorted(structure_kindname_position, key=lambda l: l[1])
        allo_kindnames_sorted = sorted(allo_kindname_position, key=lambda l: l[1])

        requires_new_kinds = not structure_kindnames_sorted == allo_kindnames_sorted
        non_magnetic = all(abs(magn) < ztol for magn in results['magnetic_moments'].get_dict().values())

        if requires_new_kinds:
            results = create_magnetic_configuration(
                structure,
                magnetic_moments,
                atol=atol,
                ztol=ztol,
                metadata={'store_provenance': False},
            )
            structure = results['results']

        return AttributeDict({
            'structure': structure,
            'magnetic_moments': None if non_magnetic else results['magnetic_moments']
        })
