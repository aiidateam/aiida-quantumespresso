# -*- coding: utf-8 -*-
from __future__ import absolute_import, division
import re
import fnmatch
import traceback

import numpy as np
from six.moves import range

from aiida.common import NotExistent, LinkType
from aiida.orm import Dict, ProjectionData, BandsData, XyData, CalcJobNode
from aiida.parsers import Parser
from aiida.plugins import OrbitalFactory

from aiida_quantumespresso.parsers import QEOutputParsingError
from aiida_quantumespresso.parsers import parse_raw_out_basic
from aiida_quantumespresso.parsers.parse_raw.projwfc import parse_lowdin_charges


def find_orbitals_from_statelines(out_file_lines, statelines_idx, positions):
    """Read all state-lines and convert into a set of orbitals.

    State-lines are lines, taken from the pseudopotential,
    describing which atomic states are used for the projection.

    :param out_file_lines: content of the projwfc stdout file
    :param statelines_idx: out_file_lines indices containing state-lines
    :param positions: list of structure site positions
    :return: orbitals, a list of orbitals suitable for setting ProjectionData
    """
    atomnum_re = re.compile(r'atom (.*?)\(')
    element_re = re.compile(r'\((.*?)\)')
    lnum_re = re.compile(r'l=(.*?)m=')
    mnum_re = re.compile(r'm=(.*?)\)')
    state_lines = [out_file_lines[wfc_line] for wfc_line in statelines_idx]
    state_dicts = []
    for state_line in state_lines:
        try:
            state_dict = {}
            state_dict['atomnum'] = int(atomnum_re.findall(state_line)[0])
            state_dict['atomnum'] -= 1  # to keep with orbital indexing
            state_dict['kind_name'] = element_re.findall(state_line)[0].strip()
            state_dict['angular_momentum'] = int(lnum_re.findall(state_line)[0])
            state_dict['magnetic_number'] = int(mnum_re.findall(state_line)[0])
            state_dict['magnetic_number'] -= 1  # to keep with orbital indexing
        except ValueError:
            raise QEOutputParsingError('State lines are not formatted in a standard way.')
        state_dicts.append(state_dict)

    # here is some logic to figure out the value of radial_nodes to use
    new_state_dicts = []
    for i in range(len(state_dicts)):  # pylint: disable=consider-using-enumerate
        radial_nodes = 0
        state_dict = state_dicts[i].copy()
        for j in range(i - 1, -1, -1):
            if state_dict == state_dicts[j]:
                radial_nodes += 1
        state_dict['radial_nodes'] = radial_nodes
        new_state_dicts.append(state_dict)
    state_dicts = new_state_dicts

    # here is some logic to assign positions based on the atom_index
    for state_dict in state_dicts:
        site_index = state_dict.pop('atomnum')
        state_dict['position'] = positions[site_index]

    # here we set the resulting state_dicts to a new set of orbitals
    orbitals = []
    real_hydrogen_orbital = OrbitalFactory('realhydrogen')
    for state_dict in state_dicts:
        orbitals.append(real_hydrogen_orbital(**state_dict))

    return orbitals


def spin_dependent_subparser(info_dict, parent_kpoints):
    """Create projection and bands arrays from the parsed data and parent calculation.

    Used to handle the different possible spin-cases in a convenient manner.

    :param info_dict: dictionary of input and parsed data
    :param parent_kpoints: kpoints from the parent calculation
    :return: ProjectionData, BandsData
    """

    out_file_lines = info_dict['out_file_lines']
    spin_down, k_states, num_bands = info_dict['spin_down'], info_dict['k_states'], info_dict['num_bands']

    #   regular expressions needed for later parsing
    regex_wavefrac1 = re.compile(r'\=(.*?)\*')  # state composition 1
    regex_wavefrac2 = re.compile(r'\+(.*?)\*')  # state comp 2
    regex_func_id = re.compile(r'\#(.*?)\]')  # state identity
    # primes arrays for the later parsing
    bands = np.zeros([k_states, num_bands])
    projection_arrays = np.zeros([k_states, num_bands, len(info_dict['wfc_lines'])])

    try:
        for i in range(k_states):
            if spin_down:
                i += k_states
            # grabs band energy
            for j in range(i * num_bands, (i + 1) * num_bands, 1):
                out_ind = info_dict['e_lines'][j]
                try:
                    # post ~6.3 output format "e ="
                    val = out_file_lines[out_ind].split()[2]
                    float(val)
                except ValueError:
                    # pre ~6.3 output format? "==== e("
                    val = out_file_lines[out_ind].split()[4]
                bands[i % k_states][j % num_bands] = val
                #subloop grabs pdos
                wave_fraction = []
                wave_id = []
                for k in range(info_dict['e_lines'][j] + 1, info_dict['psi_lines'][j], 1):
                    out_line = out_file_lines[k]
                    wave_fraction += regex_wavefrac1.findall(out_line)
                    wave_fraction += regex_wavefrac2.findall(out_line)
                    wave_id += regex_func_id.findall(out_line)
                if len(wave_id) != len(wave_fraction):
                    raise IndexError
                for l in range(len(wave_id)):  # pylint: disable=consider-using-enumerate,invalid-name
                    wave_id[l] = int(wave_id[l])
                    wave_fraction[l] = float(wave_fraction[l])
                    #sets relevant values in pdos_array
                    projection_arrays[i % k_states][j % num_bands][wave_id[l] - 1] = wave_fraction[l]
    except IndexError:
        raise QEOutputParsingError('the standard out file does not comply with the official documentation.')

    bands_data = BandsData()
    try:
        if len(info_dict['k_vect']) != len(parent_kpoints.get_kpoints()):
            raise AttributeError
        bands_data.set_kpointsdata(parent_kpoints)
    except AttributeError:
        bands_data.set_kpoints(info_dict['k_vect'].astype(float))

    bands_data.set_bands(bands, units='eV')

    orbitals = info_dict['orbitals']
    if len(orbitals) != np.shape(projection_arrays[0, 0, :])[0]:
        raise QEOutputParsingError(
            'There was an internal parsing error, '
            ' the projection array shape does not agree'
            ' with the number of orbitals'
        )
    projection_data = ProjectionData()
    projection_data.set_reference_bandsdata(bands_data)
    projections = [projection_arrays[:, :, i] for i in range(len(orbitals))]

    # Do the bands_check manually here
    for projection in projections:
        if np.shape(projection) != np.shape(bands):
            raise AttributeError('Projections not the same shape as the bands')

    # insert here some logic to assign pdos to the orbitals
    pdos_arrays = spin_dependent_pdos_subparser(info_dict)
    energy_arrays = [info_dict['energy']] * len(orbitals)
    projection_data.set_projectiondata(
        orbitals,
        list_of_projections=projections,
        list_of_energy=energy_arrays,
        list_of_pdos=pdos_arrays,
        bands_check=False
    )
    # pdos=pdos_arrays
    return bands_data, projection_data


def natural_sort_key(sort_key, _nsre=re.compile('([0-9]+)')):
    """Convert a sorting key to one suitable for natural sorting.

    For example, ``["2", "11", "1"]`` will be
    sorted to ``["1", "2", "11"]`` instead of ``["1", "11", "2"]``

    :param sort_key: Original key to be processed
    :return: list[str or int]
    """
    keys = []
    for text in _nsre.split(sort_key):
        if text.isdigit():
            keys.append(int(text))
        else:
            keys.append(text)
    return keys


def spin_dependent_pdos_subparser(out_info_dict):
    """Find and label the pdos arrays associated with the out_info_dict.

    :param out_info_dict: dictionary of input and parsed data
    :return: (pdos_name, pdos_array) tuples for all the specific pdos
    """
    spin = out_info_dict['spin']
    spin_down = out_info_dict['spin_down']
    pdos_atm_array_dict = out_info_dict['pdos_atm_array_dict']
    if spin:
        mult_factor = 2
        if spin_down:
            first_array = 4
        else:
            first_array = 3
    else:
        mult_factor = 1
        first_array = 2
    pdos_file_names = [k for k in pdos_atm_array_dict]
    pdos_file_names.sort(key=natural_sort_key)
    out_arrays = []
    # we can keep the pdos in synch with the projections by relying on the fact
    # both are produced in the same order (thus the sorted file_names)
    for name in pdos_file_names:
        this_array = pdos_atm_array_dict[name]
        for i in range(first_array, np.shape(this_array)[1], mult_factor):
            out_arrays.append(this_array[:, i])

    return out_arrays


def scan_stdout(lines):
    """Iterate through the projwf.x stdout file and record indices of required lines."""
    data = {'k_lines': [], 'e_lines': [], 'psi_lines': [], 'wfc_lines': [], 'lowdin_lines': []}
    for i, line in enumerate(lines):
        if 'k =' in line:
            data['k_lines'].append(i)
        if '==== e(' in line or line.strip().startswith('e ='):
            # The energy format in output was changed in QE6.3
            # this check supports old and new format
            data['e_lines'].append(i)
        if '|psi|^2' in line:
            data['psi_lines'].append(i)
        if 'state #' in line:
            data['wfc_lines'].append(i)
        if line.strip() == 'Lowdin Charges:':
            data['lowdin_lines'].append(i)

    return data


class ProjwfcParser(Parser):
    """
    This class is the implementation of the Parser class for projwfc.x in
    Quantum Espresso. Parses projection arrays that map the projection onto
    each point in the bands structure, as well as pdos arrays, which map
    the projected density of states onto an energy axis.
    """

    def parse(self, **kwargs):
        """
        Parses the datafolder, stores results.
        Retrieves projwfc output, and some basic information from the
        out_file, such as warnings and wall_time
        """
        # Check that the retrieved folder is there
        try:
            out_folder = self.retrieved
        except NotExistent:
            return self.exit_codes.ERROR_NO_RETRIEVED_FOLDER

        # Read standard out
        try:
            filename_stdout = self.node.get_option('output_filename')  # or get_attribute(), but this is clearer
            with out_folder.open(filename_stdout, 'r') as fil:
                out_file_lines = fil.readlines()
        except OSError:
            return self.exit_codes.ERROR_READING_OUTPUT_FILE

        job_done = False
        for line in reversed(out_file_lines):
            if 'JOB DONE' in line:
                job_done = True
                break
        if not job_done:
            self.logger.error('Computation did not finish properly')
            return self.exit_codes.ERROR_JOB_NOT_DONE

        # Parse basic info and warnings, and output them as output_parmeters
        parsed_data = parse_raw_out_basic(out_file_lines, 'PROJWFC')
        for message in parsed_data['warnings']:
            self.logger.error(message)
        output_params = Dict(dict=parsed_data)
        self.out('output_parameters', output_params)

        # check and read pdos_tot file
        out_filenames = out_folder.list_object_names()
        try:
            pdostot_filename = fnmatch.filter(out_filenames, '*pdos_tot*')[0]
            with out_folder.open(pdostot_filename, 'r') as pdostot_file:
                # Columns: Energy(eV), Ldos, Pdos
                pdostot_array = np.genfromtxt(pdostot_file)
                energy = pdostot_array[:, 0]
                dos = pdostot_array[:, 1]
        except (OSError, KeyError):
            self.logger.error('Error reading pdos_tot output file')
            return self.exit_codes.ERROR_READING_PDOSTOT_FILE

        # check and read all of the individual pdos_atm files
        pdos_atm_filenames = fnmatch.filter(out_filenames, '*pdos_atm*')
        pdos_atm_array_dict = {}
        for name in pdos_atm_filenames:
            with out_folder.open(name, 'r') as pdosatm_file:
                pdos_atm_array_dict[name] = np.genfromtxt(pdosatm_file)

        # finding the bands and projections
        # we create a dictionary the progressively accumulates more info
        out_info_dict = {'out_file_lines': out_file_lines, 'energy': energy, 'pdos_atm_array_dict': pdos_atm_array_dict}
        out_info_dict.update(scan_stdout(out_file_lines))
        try:
            new_nodes_list = self._parse_bands_and_projections(out_info_dict)
        except QEOutputParsingError as err:
            self.logger.error('Error parsing bands and projections: {}'.format(err))
            traceback.print_exc()
            return self.exit_codes.ERROR_PARSING_PROJECTIONS

        for linkname, node in new_nodes_list:
            self.out(linkname, node)

        dos_out = XyData()
        dos_out.set_x(energy, 'Energy', 'eV')
        dos_out.set_y(dos, 'Dos', 'states/eV')
        self.out('Dos', dos_out)

        try:
            node = self._parse_lodwin_charges(out_info_dict)
        except IOError as err:
            self.logger.error('Error parsing Lowdin charges: {}'.format(err))
            traceback.print_exc()
            return self.exit_codes.ERROR_PARSING_LOWDIN

        self.out('lowdin', node)

    def retrieve_parent_node(self, linkname, outgoing=True):
        """Retrieve a node from the calculation that created the parent_folder."""
        try:
            parent_calc = (
                self.node.inputs.parent_folder.get_incoming(node_class=CalcJobNode,
                                                            link_type=LinkType.CREATE).one().node
            )
        except ValueError as err:
            raise QEOutputParsingError('Could not retrieve incoming calculation of parent_folder: {}'.format(err))
        try:
            if outgoing:
                node = parent_calc.get_outgoing(link_label_filter=linkname).one().node
            else:
                node = parent_calc.get_incoming(link_label_filter=linkname).one().node
        except ValueError:
            raise QEOutputParsingError(
                'The parent calculation had no {} {} node'.format('outgoing' if outgoing else 'incoming', linkname)
            )
        return node

    def _parse_lodwin_charges(self, out_info_dict):
        data, spill_parameter = parse_lowdin_charges(out_info_dict['out_file_lines'], out_info_dict['lowdin_lines'])
        structure = self.retrieve_parent_node('structure', outgoing=False)
        try:
            site_data = [data[i + 1] for i in range(len(structure.sites))]
        except KeyError:
            raise IOError('The lowdin atom numbers do not match those in the input structure.')
        return Dict(
            dict={
                'site_data': site_data,
                'spill_parameter': spill_parameter,
                'structure_uuid': structure.uuid
            }
        )

    def _parse_bands_and_projections(self, out_info_dict):
        """Parse the standard output into bands and projection data.

        :param out_info_dict: used to pass technical internal variables to helper functions in compact form
        :return: a list of (linkname, node) containing BandsData and ProjectionData parsed from standard_out
        """
        out_file_lines = out_info_dict['out_file_lines']
        append_nodes_list = []

        # Basic check
        if len(out_info_dict['e_lines']) != len(out_info_dict['psi_lines']):
            raise QEOutputParsingError('e-lines and psi-lines are in different number')
        if len(out_info_dict['psi_lines']) % len(out_info_dict['k_lines']) != 0:
            raise QEOutputParsingError('Band Energy Points is not a multiple of kpoints')
        # calculates the number of bands
        out_info_dict['num_bands'] = len(out_info_dict['psi_lines']) // len(out_info_dict['k_lines'])

        # Use the parent input parameters, to checks if the calculation has spin.
        parent_param = self.retrieve_parent_node('output_parameters')
        spin = parent_param.get_dict().get('number_of_spin_components', 1) != 1
        out_info_dict['spin'] = spin

        # changes k-numbers to match spin
        # because if spin is on, k points double for up and down
        out_info_dict['k_states'] = len(out_info_dict['k_lines'])
        if spin:
            if out_info_dict['k_states'] % 2 != 0:
                raise QEOutputParsingError('Internal formatting error regarding spin')
            out_info_dict['k_states'] = out_info_dict['k_states'] // 2

        # adds in the k-vector for each kpoint
        k_vect = [out_file_lines[out_info_dict['k_lines'][i]].split()[2:] for i in range(out_info_dict['k_states'])]
        out_info_dict['k_vect'] = np.array(k_vect)
        structure = self.retrieve_parent_node('structure', outgoing=False)
        out_info_dict['orbitals'] = find_orbitals_from_statelines(
            out_file_lines, out_info_dict['wfc_lines'], [s.position for s in structure.sites]
        )

        kpoints = self.retrieve_parent_node('kpoints', outgoing=False)
        if spin:
            # I had to guess what the ordering of the spin is, because
            # the projwfc.x documentation doesn't say, but looking at the
            # source code I found:
            #
            # DO is=1,nspin
            #   IF (nspin==2) THEN
            #       IF (is==1) filename=trim(filproj)//'.up'
            #       IF (is==2) filename=trim(filproj)//'.down'
            #
            # Which would say that it is reasonable to assume that the
            # spin up states are written first, then spin down
            #
            out_info_dict['spin_down'] = False
            bands_data1, projection_data1 = spin_dependent_subparser(out_info_dict, kpoints)
            append_nodes_list += [('projections_up', projection_data1), ('bands_up', bands_data1)]
            out_info_dict['spin_down'] = True
            bands_data2, projection_data2 = spin_dependent_subparser(out_info_dict, kpoints)
            append_nodes_list += [('projections_down', projection_data2), ('bands_down', bands_data2)]
        else:
            out_info_dict['spin_down'] = False
            bands_data, projection_data = spin_dependent_subparser(out_info_dict, kpoints)
            append_nodes_list += [('projections', projection_data), ('bands', bands_data)]

        return append_nodes_list
