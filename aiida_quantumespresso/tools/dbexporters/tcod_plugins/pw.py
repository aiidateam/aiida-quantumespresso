# -*- coding: utf-8 -*-
"""TCOD export plugin for `PwCalculations`."""
import re

from aiida.orm.nodes.data.upf import UpfData

try:
    from aiida_tcod.tools.dbexporters.tcod import BaseTcodtranslator  # pylint: disable=import-error
except ImportError as exception:
    raise ImportError('dependency `aiida-tcod` not installed; run `pip install aiida-tcod` to do so.') from exception


class PwTcodtranslator(BaseTcodtranslator):
    """TCOD export plugin for `CpCalculations`."""

    # pylint: disable=too-many-public-methods,abstract-method

    _plugin_type_string = 'quantumespresso.pw.PwCalculation'

    _smearing_aliases = {
        'gaussian': 'Gaussian',
        'gauss': 'Gaussian',
        'methfessel-paxton': 'Methfessel-Paxton',
        'm-p': 'Methfessel-Paxton',
        'mp': 'Methfessel-Paxton',
        'marzari-vanderbilt': 'Marzari-Vanderbilt',
        'cold': 'Marzari-Vanderbilt',
        'm-v': 'Marzari-Vanderbilt',
        'mv': 'Marzari-Vanderbilt',
        'fermi-dirac': 'Marzari-Vanderbilt',
        'f-d': 'Marzari-Vanderbilt',
        'fd': 'Marzari-Vanderbilt',
    }

    _pseudopotential_types = ['NCPP', 'USPP', 'PAW']

    _upf_type_v2_regexp = re.compile(
        r"""
        \s*
        Pseudopotential type:\s*(?P<upf_type>.*)
        """, re.VERBOSE
    )

    @classmethod
    def get_software_package(cls, calc, **kwargs):  # pylint: disable=unused-argument
        """Return the package or program name that was used to produce the structure.

        Only package or program name should be used, e.g. 'VASP', 'psi3', 'Abinit', etc.
        """
        return 'Quantum ESPRESSO'

    @classmethod
    def _get_pw_energy_value(cls, calc, energy_type, **kwargs):  # pylint: disable=unused-argument
        """Return the energy of defined type in eV."""
        parameters = calc.outputs.output_parameters
        if energy_type not in parameters.attrs():
            return None
        if energy_type + '_units' not in parameters.attrs():
            raise ValueError(f'energy units for {energy_type} are unknown')
        if parameters.get_attr(energy_type + '_units') != 'eV':
            raise ValueError(
                'energy units for {} are {} instead of eV -- unit conversion is not possible yet'.format(
                    energy_type, parameters.get_attr(energy_type + '_units')
                )
            )
        return parameters.get_attr(energy_type)

    @classmethod
    def _get_atom_site_residual_force_Cartesian(cls, calc, index, **kwargs):  # pylint: disable=unused-argument,invalid-name
        """Return an array with residual force components along the Cartesian axes."""
        try:
            array = calc.outputs.output_array
            return [x[index] for x in array.get_array('forces').tolist()[-1]]
        except KeyError:
            return None

    @classmethod
    def _get_BZ_integration_grid(cls, calc, **kwargs):  # pylint: disable=unused-argument,invalid-name
        """Return an array with Brillouin zone point counts along each vector of reciprocal lattice."""
        try:
            array, _ = calc.inputs.kpoints.get_kpoints_mesh()
            return array
        except AttributeError:
            return None
        except KeyError:
            return None

    @classmethod
    def _get_BZ_integration_grid_shift(cls, calc, **kwargs):  # pylint: disable=unused-argument,invalid-name
        """Return an array with Brillouin zone point shifts along each vector of reciprocal lattice."""
        try:
            _, array = calc.inputs.kpoints.get_kpoints_mesh()
            return array
        except AttributeError:
            return None
        except KeyError:
            return None

    @classmethod
    def _get_raw_integration_smearing_method(cls, calc, **kwargs):  # pylint: disable=unused-argument
        """Return the smearing method name as string, as specified in the input parameters (if specified).

        If not 'smearing' is not
        specified, but 'occupations' == 'smearing', string with default
        value 'gaussian' is returned, as specified in
        http://www.quantum-espresso.org/wp-content/uploads/Doc/INPUT_PW.html
        """
        parameters = calc.inputs.parameters
        smearing = None
        try:
            smearing = parameters.get_dict()['SYSTEM']['smearing']
        except KeyError:
            pass
        if smearing is None:
            try:
                if parameters.get_dict()['SYSTEM']['occupations'] == 'smearing':
                    smearing = 'gaussian'
            except KeyError:
                pass
        return smearing

    @classmethod
    def get_total_energy(cls, calc, **kwargs):  # pylint: disable=unused-argument
        """Return the total energy in eV."""
        return cls._get_pw_energy_value(calc, 'energy')

    @classmethod
    def get_one_electron_energy(cls, calc, **kwargs):  # pylint: disable=unused-argument
        """Return one electron energy in eV."""
        return cls._get_pw_energy_value(calc, 'energy_one_electron')

    @classmethod
    def get_exchange_correlation_energy(cls, calc, **kwargs):  # pylint: disable=unused-argument
        """Return exchange correlation (XC) energy in eV."""
        return cls._get_pw_energy_value(calc, 'energy_xc')

    @classmethod
    def get_ewald_energy(cls, calc, **kwargs):  # pylint: disable=unused-argument
        """Return Ewald energy in eV."""
        return cls._get_pw_energy_value(calc, 'energy_ewald')

    @classmethod
    def get_hartree_energy(cls, calc, **kwargs):  # pylint: disable=unused-argument
        """Return Hartree energy in eV."""
        return cls._get_pw_energy_value(calc, 'energy_hartree')

    @classmethod
    def get_fermi_energy(cls, calc, **kwargs):  # pylint: disable=unused-argument
        """Return Fermi energy in eV."""
        return cls._get_pw_energy_value(calc, 'fermi_energy')

    @classmethod
    def get_number_of_electrons(cls, calc, **kwargs):  # pylint: disable=unused-argument
        """Return the number of electrons."""
        parameters = calc.outputs.output_parameters
        if 'number_of_electrons' not in parameters.attrs():
            return None
        return parameters.get_attr('number_of_electrons')

    @classmethod
    def get_computation_wallclock_time(cls, calc, **kwargs):  # pylint: disable=unused-argument
        """Return the computation wallclock time in seconds."""
        parameters = calc.outputs.output_parameters
        if 'wall_time_seconds' not in parameters.attrs():
            return None
        return parameters.get_attr('wall_time_seconds')

    @classmethod
    def get_atom_site_residual_force_Cartesian_x(cls, calc, **kwargs):  # pylint: disable=unused-argument,invalid-name
        """Return a list of x components for Cartesian coordinates of residual force for atom.

        The list order MUST be the same as in the resulting structure.
        """
        return cls._get_atom_site_residual_force_Cartesian(calc, 0)

    @classmethod
    def get_atom_site_residual_force_Cartesian_y(cls, calc, **kwargs):  # pylint: disable=unused-argument,invalid-name
        """Return a list of y components for Cartesian coordinates of residual force for atom.

        The list order MUST be the same as in the resulting structure.
        """
        return cls._get_atom_site_residual_force_Cartesian(calc, 1)

    @classmethod
    def get_atom_site_residual_force_Cartesian_z(cls, calc, **kwargs):  # pylint: disable=unused-argument,invalid-name
        """Return a list of z components for Cartesian coordinates of residual force for atom.

        The list order MUST be the same as in the resulting structure.
        """
        return cls._get_atom_site_residual_force_Cartesian(calc, 2)

    @classmethod
    def get_BZ_integration_grid_X(cls, calc, **kwargs):  # pylint: disable=invalid-name
        """Return a number of points in the Brillouin zone along reciprocal lattice vector X."""
        array = cls._get_BZ_integration_grid(calc, **kwargs)
        if array is not None:
            return array[0]
        return None

    @classmethod
    def get_BZ_integration_grid_Y(cls, calc, **kwargs):  # pylint: disable=invalid-name
        """Return a number of points in the Brillouin zone along reciprocal lattice vector Y."""
        array = cls._get_BZ_integration_grid(calc, **kwargs)
        if array is not None:
            return array[1]

    @classmethod
    def get_BZ_integration_grid_Z(cls, calc, **kwargs):  # pylint: disable=invalid-name
        """Return a number of points in the Brillouin zone along reciprocal lattice vector Z."""
        array = cls._get_BZ_integration_grid(calc, **kwargs)
        if array is not None:
            return array[2]
        return None

    @classmethod
    def get_BZ_integration_grid_shift_X(cls, calc, **kwargs):  # pylint: disable=invalid-name
        """Return the shift of the Brillouin zone points along reciprocal lattice vector X."""
        array = cls._get_BZ_integration_grid_shift(calc, **kwargs)
        if array is not None:
            return array[0]

    @classmethod
    def get_BZ_integration_grid_shift_Y(cls, calc, **kwargs):  # pylint: disable=invalid-name
        """Return the shift of the Brillouin zone points along reciprocal lattice vector Y."""
        array = cls._get_BZ_integration_grid_shift(calc, **kwargs)
        if array is not None:
            return array[1]

    @classmethod
    def get_BZ_integration_grid_shift_Z(cls, calc, **kwargs):  # pylint: disable=invalid-name
        """Return the shift of the Brillouin zone points along reciprocal lattice vector Z."""
        array = cls._get_BZ_integration_grid_shift(calc, **kwargs)
        if array is not None:
            return array[2]

    @classmethod
    def get_integration_smearing_method(cls, calc, **kwargs):
        """Return the smearing method name as string."""
        smearing = cls._get_raw_integration_smearing_method(calc, **kwargs)
        if smearing is None:
            return None
        if smearing in cls._smearing_aliases:
            return cls._smearing_aliases[smearing]
        return 'other'

    @classmethod
    def get_integration_smearing_method_other(cls, calc, **kwargs):  # pylint: disable=invalid-name
        """Return the smearing method name as string if the name is different from specified in cif_dft.dic."""
        smearing = cls._get_raw_integration_smearing_method(calc, **kwargs)
        if smearing is None or smearing in cls._smearing_aliases:
            return None
        return smearing

    @classmethod
    def get_integration_Methfessel_Paxton_order(cls, calc, **kwargs):  # pylint: disable=invalid-name
        """Return the order of Methfessel-Paxton approximation if used."""
        if cls.get_integration_smearing_method(calc, **kwargs) == 'Methfessel-Paxton':
            return 1

    @classmethod
    def get_kinetic_energy_cutoff_wavefunctions(cls, calc, **kwargs):  # pylint: disable=invalid-name,unused-argument
        """Return kinetic energy cutoff for wavefunctions in eV."""
        from qe_tools import CONSTANTS
        parameters = calc.inputs.parameters
        ecutwfc = None
        try:
            ecutwfc = parameters.get_dict()['SYSTEM']['ecutwfc']
        except KeyError:
            pass
        if ecutwfc is None:
            return None
        return ecutwfc * CONSTANTS.ry_to_ev

    @classmethod
    def get_kinetic_energy_cutoff_charge_density(cls, calc, **kwargs):  # pylint: disable=invalid-name,unused-argument
        """Return kinetic energy cutoff for charge density in eV.

        .. note :: by default returns 4 * ecutwfc, as indicated in
            http://www.quantum-espresso.org/wp-content/uploads/Doc/INPUT_PW.html
        """
        from qe_tools import CONSTANTS
        parameters = calc.inputs.parameters
        try:
            return parameters.get_dict()['SYSTEM']['ecutrho'] * CONSTANTS.ry_to_ev
        except KeyError:
            pass
        ecutwfc = cls.get_kinetic_energy_cutoff_wavefunctions(calc)
        if ecutwfc is None:
            return None
        return 4 * ecutwfc

    @classmethod
    def get_kinetic_energy_cutoff_EEX(cls, calc, **kwargs):  # pylint: disable=invalid-name,unused-argument
        """Return kinetic energy cutoff for exact exchange (EEX) operator in eV.

        .. note :: by default returns ecutrho, as indicated in
            http://www.quantum-espresso.org/wp-content/uploads/Doc/INPUT_PW.html
        """
        from qe_tools import CONSTANTS
        parameters = calc.inputs.parameters
        try:
            return parameters.get_dict()['SYSTEM']['ecutfock'] * CONSTANTS.ry_to_ev
        except KeyError:
            pass
        return cls.get_kinetic_energy_cutoff_charge_density(calc)

    @classmethod
    def _get_raw_pseudopotential_type(cls, calc, **kwargs):  # pylint: disable=unused-argument
        """Return the raw pseudopotential type."""
        from aiida.orm.nodes.data.upf import parse_upf
        types = {}
        for node in calc.get_incoming(node_class=UpfData).all_nodes():
            element = node.element
            parsed = parse_upf(node.get_file_abs_path())
            if parsed['version'] != '2':
                types[element] = None
                continue
            upf_type = None
            with open(node.filename) as handle:
                for line in handle:
                    match = cls._upf_type_v2_regexp.match(line.strip())
                    if match:
                        upf_type = match.group('upf_type')
                        break
            types[element] = upf_type
        return types

    @classmethod
    def get_pseudopotential_atom_type(cls, calc, **kwargs):
        """Return a list of atom types.

        Each atom type MUST occur only once in this list. List MUST be sorted.
        """
        raw_types = cls._get_raw_pseudopotential_type(calc, **kwargs)
        return sorted(raw_types)

    @classmethod
    def get_pseudopotential_type(cls, calc, **kwargs):
        """Return a list of pseudopotential types.

        List MUST be sorted by atom types.
        """
        types = []
        raw_types = cls._get_raw_pseudopotential_type(calc, **kwargs)
        for element in sorted(raw_types):
            if raw_types[element] is None or raw_types[element] in cls._pseudopotential_types:
                types.append(raw_types[element])
            else:
                types.append('other')
        return types

    @classmethod
    def get_pseudopotential_type_other_name(cls, calc, **kwargs):  # pylint: disable=invalid-name
        """Return a list of other pseudopotential type names.

        List MUST be sorted by atom types.
        """
        types = []
        raw_types = cls._get_raw_pseudopotential_type(calc, **kwargs)
        for element in sorted(raw_types):
            if raw_types[element] is None or raw_types[element] in cls._pseudopotential_types:
                types.append(None)
            else:
                types.append(raw_types[element])
        return types
